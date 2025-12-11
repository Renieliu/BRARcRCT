library(shiny)
library(rstan)
library(MASS)
library(tidyr)
library(ggplot2)
library(dplyr)
library(shinyBS)
library(rsconnect)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ui <- fluidPage(
  titlePanel("BRAR cRCT"),
  
# Description paragraph below the title
p("This app simulates Bayesian Response-Adaptive Randomization (BRAR) in a cluster randomized controlled trial (cRCT). 
You can customize trial setup parameters, run simulations, and obtain suggested treatment assignments for each cluster. 
Cluster sizes can vary across clusters, but you must specify the maximum possible cluster size at the beginning. 
The outcome is assumed to be continuous, and higher values are interpreted as better treatment effects."),
  
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "output.phase == 'setup'",
        numericInput("n_cluster", "Number of clusters in total:", 50, min = 2),
        bsTooltip("n_cluster", 
                  "Total number of clusters enrolled in the trial. This determines the sample size at the cluster level.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
        
        numericInput("n_each_max", "Maximum number of subjects per cluster:", 8, min = 2),
        bsTooltip("n_each_max", 
                  "Maximum number of individual subjects allowed within each cluster.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
        
        numericInput("n_repeats", "Repeats per cluster in the simulation:", 500, min = 10),
        bsTooltip("n_repeats", 
                  "Number of repeated simulations for each cluster to obtain the assignment distribution.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
        
        numericInput("initial", "Number of initial clusters:", 4, min = 1),
        bsTooltip("initial", 
                  "Number of clusters initially enrolled before adaptive allocation begins.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
  
        numericInput("U", "Upper threshold (U):", 0.98),
        bsTooltip("U", 
                  "Posterior probability threshold for decision making.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
        
        numericInput("L", "Lower threshold (L):", 0.02),
        bsTooltip("L", 
                  "Posterior probability threshold for decision making.",
                  placement = "right", trigger = "hover", options = list(container = "body")),
        
        actionButton("start", "Start Simulation")
      ),
      
      conditionalPanel(
        condition = "output.phase != 'setup'",
        uiOutput("phase"),
        uiOutput("inputs_ui")
      )
    ),
    
    mainPanel(
      plotOutput("allocation_plot"),
      textOutput("next_assignment"),
      textOutput("decision_text"),
      textOutput("final_assignments")
    )
  )
)

server <- function(input, output, session) {
  
  # initialization
  
  rv <- reactiveValues(
    n_cluster = NULL, n_each = NULL, n_repeats = NULL, initial = NULL,
    y_A = list(),
    y_B = list(),
    outcomes = list(), allocations = NULL,
    next_cluster = NULL,
    plot_data = NULL,
    next_allocation_dist = NULL,
    next_assigned = NULL,
    decision_text = NULL,
    setup_done = FALSE,
    current_cluster_size = NULL
  )
  
  # 1. render phase text (optional UI display only)
  
  output$phase <- renderText({
    if (!rv$setup_done) {
      "setup"
    } else if (rv$next_cluster <= rv$initial) {
      paste("Inital stage - Cluster", rv$next_cluster)
    } else if (rv$next_cluster <= rv$n_cluster) {
      paste("Adaptive stage - Cluster", rv$next_cluster, "(suggested:", rv$next_assigned, ")")
    } else {
      "All clusters done!"
    }
  })
  outputOptions(output, "phase", suspendWhenHidden = FALSE)
  
  # 2. initialize parameters and reset state
  
  observeEvent(input$start, {
    rv$n_cluster <- input$n_cluster
    rv$n_each_max <- input$n_each_max
    rv$n_repeats <- input$n_repeats
    rv$initial <- input$initial
    
    rv$outcomes <- vector("list", rv$n_cluster)
    rv$cluster_sizes <- rep(NA, rv$n_cluster)
    rv$y_A <- list()
    rv$y_B <- list()
    rv$allocations <- rep(NA, rv$n_cluster)
    rv$next_cluster <- 1
    rv$plot_data <- NULL
    rv$next_allocation_dist <- NULL
    rv$next_assigned <- "A"
    rv$decision_text <- NULL
    rv$setup_done <- TRUE
    rv$current_cluster_size <- NULL
  })
  
  # 3. observe confirm_cluster_size here
  
  observeEvent(input$confirm_cluster_size, {
    if (rv$setup_done && rv$next_cluster <= rv$n_cluster) {
      rv$current_cluster_size <- input$n_each_current
      rv$cluster_sizes[rv$next_cluster] <- input$n_each_current
    }
  })
  
  # 4. render UI inputs dynamically based on cluster size
  
  output$inputs_ui <- renderUI({
    if (rv$setup_done && rv$next_cluster <= rv$n_cluster) {
      if (is.null(rv$current_cluster_size)) {
        tagList(
          numericInput("n_each_current",
                       paste("Enter number of subjects in cluster", rv$next_cluster),
                       value = rv$n_each_max,
                       min = 1, max = rv$n_each_max),
          actionButton("confirm_cluster_size", "Confirm Cluster Size")
        )
        
      } else {
        inputs <- lapply(1:rv$current_cluster_size, function(i) {
          numericInput(paste0("y_", i), paste("Outcome for participant", i), value = NA)
        })
        
        inputs <- append(inputs, list(
          radioButtons("group_choice",
                       if (rv$next_cluster <= rv$initial) {
                         "Select group for this cluster:"
                       } else {
                         paste0("Suggested: ", rv$next_assigned, ". Choose final group:")
                       },
                       choices = c("A", "B"),
                       selected = rv$next_assigned,
                       inline = TRUE),
          actionButton("submit", "Submit Outcomes")
        ))
        
        return(inputs)
      }
    }
  })
  
  # 5. process submitted outcomes
  
  observeEvent(input$submit, {
    y_input <- matrix(sapply(1:rv$current_cluster_size, function(i) input[[paste0("y_", i)]]), nrow = 1)
    
    
    assigned <- input$group_choice
    rv$allocations[rv$next_cluster] <- assigned
    
    rv$outcomes[[rv$next_cluster]] <- y_input
    if (assigned == "A") {
      rv$y_A[[length(rv$y_A) + 1]] <- y_input
    } else {
      rv$y_B[[length(rv$y_B) + 1]] <- y_input
    }
    
    rv$current_cluster_size <- NULL
    
    if (rv$next_cluster < rv$initial) {
      ## Burn-in → just advance
      rv$next_cluster <- rv$next_cluster + 1
      
      if (rv$next_cluster == rv$initial + 1) {
#        withProgress(message = "Running BRAR adaptive simulation...", value = 0, {
#          incProgress(0.1, detail = "Preparing data...")
          
          stan_data <- list(
            N_A = length(rv$y_A),
            N_B = length(rv$y_B),
            cluster_size_A = sapply(rv$y_A, length),
            cluster_size_B = sapply(rv$y_B, length),
            cluster_size_max = rv$n_each_max,
            Y_A = lapply(rv$y_A, function(x) c(x, rep(0, rv$n_each_max - length(x)))),
            Y_B = lapply(rv$y_B, function(x) c(x, rep(0, rv$n_each_max - length(x)))),
            A_sigma_W = 1,  # Half-Cauchy scale based on expected var
            rho_alpha = 1,                 # Beta(1,1) = Uniform(0,1) on ICC
            rho_beta  = 1 
          )
          
#          incProgress(0.3, detail = "Sampling from posterior...")
          fit <- stan(file = "BRAR_cRCT_vary.stan", data = stan_data,
                      iter = 2000, chains = 4, refresh = 100)
          
#          incProgress(0.3, detail = "Extracting results...")
          post <- rstan::extract(fit)
          
#          incProgress(0.2, detail = "Simulating allocation...")
          allocation_dist <- replicate(rv$n_repeats, {
            mu_A <- sample(post$mu_A, 1)
            mu_B <- sample(post$mu_B, 1)
            ifelse(mu_A > mu_B, "A", "B")
          })
          
          rv$next_allocation_dist <- allocation_dist
          rv$next_assigned <- sample(allocation_dist, 1)
          
          df <- data.frame(Cluster = rv$next_cluster, Group = allocation_dist)
          rv$plot_data <- bind_rows(rv$plot_data, df)
          
          prop_B <- mean(allocation_dist == "B")
          if (prop_B > input$U) {
            rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Stop for efficacy", prop_B)
          } else if (prop_B < input$L) {
            rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Stop for futility", prop_B)
          } else {
            rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Continue", prop_B)
          }
#        })
      }
      
    } else if (rv$next_cluster < rv$n_cluster) {
#      withProgress(message = "Running BRAR adaptive simulation...", value = 0, {
#        incProgress(0.1, detail = "Preparing data...")
        
        stan_data <- list(
          N_A = length(rv$y_A),
          N_B = length(rv$y_B),
          cluster_size_A = sapply(rv$y_A, length),
          cluster_size_B = sapply(rv$y_B, length),
          cluster_size_max = rv$n_each_max,
          Y_A = lapply(rv$y_A, function(x) c(x, rep(0, rv$n_each_max - length(x)))),
          Y_B = lapply(rv$y_B, function(x) c(x, rep(0, rv$n_each_max - length(x)))),
          A_sigma_W = 1,  # Half-Cauchy scale based on expected var
          rho_alpha = 1,                 # Beta(1,1) = Uniform(0,1) on ICC
          rho_beta  = 1 
        )
        
#        incProgress(0.3, detail = "Sampling from posterior...")
        fit <- stan(file = "BRAR_cRCT_vary.stan", data = stan_data,
                    iter = 2000, chains = 4, refresh = 100)
        
#        incProgress(0.3, detail = "Extracting results...")
        post <- rstan::extract(fit)
        
#        incProgress(0.2, detail = "Simulating allocation...")
        allocation_dist <- replicate(rv$n_repeats, {
          mu_A <- sample(post$mu_A, 1)
          mu_B <- sample(post$mu_B, 1)
          ifelse(mu_A > mu_B, "A", "B")
        })
        
        rv$next_allocation_dist <- allocation_dist
        rv$next_assigned <- sample(allocation_dist, 1)
        
        df <- data.frame(Cluster = rv$next_cluster, Group = allocation_dist)
        rv$plot_data <- bind_rows(rv$plot_data, df)
        
        prop_B <- mean(allocation_dist == "B")
        if (prop_B > input$U) {
          rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Stop for efficacy", prop_B)
        } else if (prop_B < input$L) {
          rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Stop for futility", prop_B)
        } else {
          rv$decision_text <- sprintf("Posterior probability of treatment B better than treatment A = %.3f → Continue", prop_B)
        }
#      })
      rv$next_cluster <- rv$next_cluster + 1
      #      rv$current_cluster_size <- NULL
    } else {
      ## Final cluster → done
      rv$next_cluster <- rv$next_cluster + 1
    }
  })
  
  # 6. render plot + text outputs
  
  output$allocation_plot <- renderPlot({
    if (is.null(rv$plot_data)) return(NULL)
    
    df_summary <- rv$plot_data %>%
      group_by(Cluster, Group) %>%
      summarise(Count = n(), .groups = "drop") %>%
      mutate(Proportion = Count / rv$n_repeats)
    
    ggplot(df_summary, aes(x = factor(Cluster), y = Proportion, fill = Group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("A" = "#FDBE85", "B" = "#6BAED6")) +
      labs(x = "Cluster", y = "Proportion", fill = "Group") +
      theme_minimal()
  })
  
  output$next_assignment <- renderText({
    if (!is.null(rv$next_assigned) && rv$next_cluster > rv$initial && rv$next_cluster <= rv$n_cluster) {
      paste("Suggested group for cluster", rv$next_cluster, ":", rv$next_assigned)
    }
  })
  
  output$decision_text <- renderText({
    rv$decision_text
  })
  
  output$final_assignments <- renderText({
    if (rv$setup_done && rv$next_cluster > rv$n_cluster) {
      prop_A <- mean(rv$allocations[-(1:rv$initial)] == "A", na.rm = TRUE)
      prop_B <- 1 - prop_A
      sprintf("Final assignments (excluding initial): %.1f%% A, %.1f%% B",
              prop_A * 100, prop_B * 100)
    }
  })
  
}

shinyApp(ui, server)
