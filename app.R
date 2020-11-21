# Load packages
suppressPackageStartupMessages({
    library(Biostrings)
    library(magrittr)
    library(shiny)
    library(shinythemes)
    library(tidyverse)
    library(extrafont)
})

# Load data
gene <- DNAString("ATGTGGGCTCTCCTGAGAAGGCATTGGCAGTCAAACAGTTACAGTTGGCAGGAAGGGTTCCAGGCATTTCTCCACCAGATGAAAATGGCGTAAGGACATATCATAGACCCGAGGCCAAGTCAGGAGGTCATGGTGTTGGCTGGAGTCCAATGATTCCATACTCATTTAAAGTACCTCCAGAATGGGAAGAGGTTCCTGTTTCAATTGCAGATCTTGGTGGCACAGAGATTGATTTGAGGTTTGTAAACCAAAAGGAAGGAAGCTTGTCTGTTGTTGTTGCACCGGTTCTAAGATTTGCAGACAATGATGAAAACGCTAGAATTGAAGACATTGGGCCACCTGAGAAAGTGATATATGCTTTTGGCCCAGAGCTGATTGGAGATAATATAGAAGGGAAGGTAAAAAGTGCTGAAGTCAAGGAGCAATCTGGTAGAAAATATTATCAATATGAGATAGAATCACCACATGCACTAGTCACTGCGACGGCTGCTGGCAATAGGATATATTTGATGACTGTCAGTGCAAGTGGTCTTCAATGGAAAAAACATTTTCCAAATTTAAAGAAGATTGCGGACTCCTTCCGAGTGGATTAG")
cuts <- matchPattern("GG",gene)
df <- tibble(StartPosition=start(cuts),
             EndPosition=end(cuts),
             NucleicSequence=as.character(resize(cuts,3,"end")))
hpal <- colorRampPalette(c("red","white","green"))(length(cuts))
ChristmasTree <- read.csv("https://raw.githubusercontent.com/t-redactyl/Blog-posts/master/Christmas%20tree%20base%20data.csv")

Desired.Lights <- 50
Total.Lights <- sum(round(Desired.Lights * 0.35) + round(Desired.Lights * 0.20) + 
                        round(Desired.Lights * 0.17) + round(Desired.Lights * 0.13) +
                        round(Desired.Lights * 0.10) + round(Desired.Lights * 0.05))

Lights <- data.frame(Lights.X = c(round(runif(round(Desired.Lights * 0.35), 4, 18), 0),
                                  round(runif(round(Desired.Lights * 0.20), 5, 17), 0),
                                  round(runif(round(Desired.Lights * 0.17), 6, 16), 0),
                                  round(runif(round(Desired.Lights * 0.13), 7, 15), 0),
                                  round(runif(round(Desired.Lights * 0.10), 8, 14), 0),
                                  round(runif(round(Desired.Lights * 0.05), 10, 12), 0)))
Lights$Lights.Y <- c(round(runif(round(Desired.Lights * 0.35), 4, 6), 0),
                     round(runif(round(Desired.Lights * 0.20), 7, 8), 0),
                     round(runif(round(Desired.Lights * 0.17), 9, 10), 0),
                     round(runif(round(Desired.Lights * 0.13), 11, 12), 0),
                     round(runif(round(Desired.Lights * 0.10), 13, 14), 0),
                     round(runif(round(Desired.Lights * 0.05), 15, 17), 0))
Lights$Lights.Colour <- c(round(runif(Total.Lights, 1, 4), 0))

Baubles <- data.frame(Bauble.X = c(6, 9, 15, 17, 5, 13, 16, 7, 10, 14, 7, 9, 11, 
                                   14, 8, 14, 9, 12, 11, 12, 14, 11, 17, 10))
Baubles$Bauble.Y <- c(4, 5, 4, 4, 5, 5, 5, 6, 6, 6, 8, 8, 8, 8, 10,
                      10, 11, 11, 12, 13, 10, 16, 7, 14)
Baubles$Bauble.Colour <- factor(c(1, 2, 2, 3, 2, 3, 1, 3, 1, 1, 1, 2, 1, 2,
                                  3, 3, 2, 1, 3, 2, 1, 3, 3, 1))
Baubles$Bauble.Size <- c(1, 3, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 3, 3, 3,
                         2, 3, 1, 1, 2, 2, 3, 3, 2)

#font_import()
loadfonts(quiet = TRUE)

# plot function
xmasTree <- function(col=hpal[length(cuts)]){
    # from https://t-redactyl.io/blog/2016/12/a-very-ggplot2-christmas.html
    
    ChristmasTree %<>% mutate(
        Tree.Colour=replace(Tree.Colour,Tree.Y > 3, col))
    
    tree <- ggplot() + 
        geom_tile(data = ChristmasTree, aes(x = Tree.X, 
                                            y = Tree.Y, 
                                            fill = Tree.Colour)) +       
        scale_fill_identity() + 
        theme_bw() +
        scale_x_continuous(breaks = NULL) + 
        scale_y_continuous(breaks = NULL) +
        labs(x = "", y = "")
    
    tree <- tree +
        geom_point(data = Lights, aes(x = Lights.X, y = Lights.Y, alpha = Lights.Colour),
                   colour = "lightgoldenrodyellow", shape = 16) +
        theme(legend.position = "none")
 
    tree <- tree + 
        geom_point(data = Baubles, aes(x = Bauble.X, y = Bauble.Y, 
                                       colour = Bauble.Colour, size = Bauble.Size),
                   shape = 16) +
        scale_colour_manual(values = c("firebrick2", "gold", "dodgerblue3")) +
        scale_size_area(max_size = 12)
    
    tree <- tree +
        geom_segment(aes(x = 2.5, xend = 4.5, y = 1.5, yend = 1.5), colour = "blueviolet", size = 2) +
        geom_segment(aes(x = 5.5, xend = 8.5, y = 1.5, yend = 1.5), colour = "dodgerblue3", size = 2) +
        geom_segment(aes(x = 13.5, xend = 16.5, y = 1.5, yend = 1.5), colour = "blueviolet", size = 2) +
        geom_segment(aes(x = 17.5, xend = 19.5, y = 1.5, yend = 1.5), colour = "dodgerblue3", size = 2) +
        geom_segment(aes(x = 3.5, xend = 3.5, y = 0.5, yend = 2.5), colour = "blueviolet", size = 2) +
        geom_segment(aes(x = 7.0, xend = 7.0, y = 0.5, yend = 2.5), colour = "dodgerblue3", size = 2) +
        geom_segment(aes(x = 15.0, xend = 15.0, y = 0.5, yend = 2.5), colour = "blueviolet", size = 2) +
        geom_segment(aes(x = 18.5, xend = 18.5, y = 0.5, yend = 2.5), colour = "dodgerblue3", size = 2)
    
    tree <- tree +
        annotate("text", x = 11, y = 20, label = "Merry Christmas!", 
                 family = "Luminari", size = 12)
    return(tree)
}

# colorname function
colName <- function(col=hpal[length(cuts)]){
    # from https://gist.github.com/sklarz-bgu/01a550f59cdf5bc85a48e15f5e94a6ba
    colors()                 %>%        # Get X11 colors (hex)
        col2rgb              %>%        # Convert to RGB matrix
        data.frame           %>%        # Convert to data.frame
        setNames(.,colors()) %>%        # Set color names
        t                    %>%        # Transform so colors are in rows (columns: R,G,B)
        data.frame           %>%        # Re-convert to data.frame
        apply(.,1,function(x) sum(abs(x-(col %>% col2rgb()))) ) %>%  # For each color, calc the sum of diff between mycol RGB and the color RGB
        sort                 %>%        # Sort so color with smallest diff comes first
        '['(1)               %>%        # Get the first color in the df = closest
        names %>% return()
}

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                titlePanel("Val din julgran färg!"),
                div(HTML(paste("Select one of the position in the slider below",
                          "<br/>that can be cut with the gene scissors.",
                          "<br/>Observe the corresponding color change!"))),
                sidebarLayout(
                    sidebarPanel(
                        # Select one of the cutting site
                        sliderInput(inputId = "cutSite",
                                    label="Select a cut site",
                                    min = 1,
                                    max = length(cuts),
                                    step = 1,
                                    value = length(cuts))
                    ),
                    # Show a plot of the generated distribution
                    mainPanel(
                        htmlOutput(outputId = "tableDesc"),
                        tableOutput("tab"),
                        htmlOutput(outputId = "plotDesc"),
                        plotOutput("treePlot"),
                        htmlOutput(outputId = "Disclaimer")
                    )
                )
)

# Define server function
server <- function(input, output) {
    
    output$tableDesc <- renderUI(HTML(
        paste("You have selected the cut site number",
              "<b>",input$cutSite,"</b>",
              "<br/>In the table below you can see where it starts and ends in the gene, as well as the recognised sequence (guide RNA)")))
    output$tab <- renderTable(df[input$cutSite,],
                              bordered = TRUE,align = "c")
    # get the color name
    output$plotDesc <- renderUI(HTML(
        paste("The truncated protein now encodes for a",
              "<b>",colName(col=hpal[input$cutSite]),"</b>","pigment")))
    
    output$treePlot <- renderPlot({
        xmasTree(col=hpal[input$cutSite])
    })
    
    output$Disclaimer <- renderUI(
        HTML(paste(
            "<h3>Disclaimer</h3>",
            "<br/>There is no such known protein that would change the tree color as described",
                 "<br/>While this example is purely <b>fictional</b>",
                 "the technology described exist and could be used as described.",
                 "<br/><br/>Some code used for this webpage was obtained from",
                 "<a href='https://t-redactyl.io/blog/2016/12/a-very-ggplot2-christmas.html'>here</a>",
                 "and <a href='https://gist.github.com/sklarz-bgu/01a550f59cdf5bc85a48e15f5e94a6ba'>there</a>."
    )))
}
# Create Shiny object
shinyApp(ui = ui, server = server)
