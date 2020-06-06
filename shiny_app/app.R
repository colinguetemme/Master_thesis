library(shiny)
library(shinyMatrix)
library(expm)
library(phasty)
library(shinydashboard)
library(ggplot2)

sub <- matrix(0.2, 3, 3)
init <- matrix(c(1, 0, 0), 1, 3)
reward <- matrix(c(1, 1, 1), 1, 3)

ui <- dashboardPage(
    dashboardHeader(title = 'Phasty app'),
    dashboardSidebar(
        sidebarMenu(
            menuItem('Plot univariate', tabName = 'univ'),
            menuItem('Reward', tabName = 'reward')
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = 'univ',
                    splitLayout(
                        tabsetPanel(
                            tabPanel('pdf', plotOutput("plot1")),
                            tabPanel('cdf', plotOutput('plot2')),
                            tabPanel('quantile', plotOutput('plot3')),
                            tabPanel('summary', verbatimTextOutput('test')),
                            id = ('menu_plot')
                        ),
                        verticalLayout(
                            fluidRow(
                                column(2, actionButton('add', 'add'), offset = 1),
                                column(2, actionButton('remove', 'remove')),
                                column(2, actionButton('clear', 'clear'), offset = 1)
                            ),
                            fluidRow(column(8, h2('subintensity matrix'), offset = 1)),
                            fluidRow(
                                column(8, matrixInput('subint_mat', value <- sub), offset = 1)
                            ),
                            fluidRow(column(8, h2('initial probabilities'), offset = 1)),
                            fluidRow(
                                column(8, matrixInput('init_probs', value <- init), offset = 1)
                            ),
                            fluidRow(column(6, h2('Rewards'), offset = 1)),
                            fluidRow(
                                column(6, matrixInput('reward', value <- reward), offset = 1),
                                column(1, actionButton('apply', 'apply'), offset = 1)
                            )
                        )
                        
                    )
            ),
            tabItem(tabName = 'reward',
                    fluidRow(column(10, h2('subintensity matrix'), offset = 1)),
                    fluidRow(
                        column(10, matrixInput('subint_mat', value <- sub), offset = 1)
                    ),
                    fluidRow(column(10, h2('initial probabilities'), offset = 1)),
                    fluidRow(
                        column(10, matrixInput('init_probs', value <- init), offset = 1)
                    ),
                    fluidRow(column(10, h2('initial probabilities'), offset = 1)),
                    fluidRow(
                        column(10, matrixInput('reward', value <- init), offset = 1)
                    )
            )
        )
        
    )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output, session) {
    
    data <- reactive(phase_type(input$subint_mat, input$init_probs, 6))
    
    max_x <- reactive(qphtype(0.95, data()))
    
    pas <- reactive({
        if(class(data()) == 'disc_phase_type'){
            pas <- 1
        } else if (class(data()) == 'cont_phase_type'){
            pas = 0.01 * max_x()
        }
    })
    
    x <- reactive({
        if (input$menu_plot == 'pdf'){
            x = seq(0, max_x(), pas())
        } else if (input$menu_plot == 'cdf'){
            x = seq(0, max_x(), pas())
        } else if (input$menu_plot == 'quantile'){
            x = seq(0.03, 0.97, 0.01)
        }
    })
    
    output$plot1 <- renderPlot({
        ploty <- cbind(as.matrix(x()), as.matrix(dphtype(x(), data())))
        ploty <-  data.frame(ploty)
        ggplot(ploty) + geom_area(aes(x = X1, y = X2), fill = 'red', alpha = 0.6)
    })
    output$plot2 <- renderPlot({
        ploty <- cbind(as.matrix(x()), as.matrix(pphtype(x(), data())))
        ploty <-  data.frame(ploty)
        ggplot(ploty) + geom_area(aes(x = X1, y = X2), fill = 'red', alpha = 0.6)
    })
    output$plot3 <- renderPlot({
        ploty <- cbind(as.matrix(x()), as.matrix(qphtype(x(), data())))
        ploty <-  data.frame(ploty)
        ggplot(ploty) + geom_area(aes(x = X1, y = X2), fill = 'red', alpha = 0.6)
    })
    
    output$test <- renderPrint(summary(data()))
    
    observeEvent(input$add, {
        updateMatrixInput(session, "subint_mat",
                          value = rbind(cbind(input$subint_mat, 0), 0))
        updateMatrixInput(session, "init_probs",
                          value = cbind(input$init_probs, 0))
    })
    observeEvent(input$remove, {
        updateMatrixInput(session, "subint_mat",
                          value = input$subint_mat[1:(nrow(input$subint_mat)-1),
                                                   -ncol(input$subint_mat)])
        updateMatrixInput(session, "init_probs", value = t(matrix(
            input$init_probs[1, -ncol(input$init_probs)])))
    })
    observeEvent(input$clear, {
        updateMatrixInput(session, "subint_mat", value = 0 * matrix(
            as.numeric(input$subint_mat), ncol = ncol(input$subint_mat)))
        updateMatrixInput(session, 'init_probs',
                          value = t(matrix(c(1, rep(0, ncol(input$init_probs)-1)))))
    })
}

shinyApp(ui, server)
