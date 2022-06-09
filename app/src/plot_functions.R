plot_time_series <- function(input, klass, koos, plot_nr) {

      
    # Change linecolor in case of dark mode
    if (input$dark) {
      line_color="grey80"
    } else if (input$dark==FALSE) {
      line_color="black"
    } 
    
    ############## Adjust X-scaling if necessary ##########
    
    #Adjust scale if range for x (min,max) is specified
    if (plot_nr == 1) {
    
        if (input$range_x != "" &&  input$change_scale == TRUE && !vals$Datum) {
            rng_x <- as.numeric(strsplit(input$range_x,",")[[1]])
            observe({ print(rng_x) })
      
      
            #If min>max invert the axis
            if (rng_x[1]>rng_x[2]) {
                p <- p+ scale_x_reverse()
                klaas <-  klaas %>% filter(Time >= rng_x[2] & Time <= rng_x[1] )
                koos <-  koos %>% filter(Time >= rng_x[2] & Time <= rng_x[1] )
            } else {
            #Select timepoints within the set limits
                klaas <-  klaas %>% filter(Time >= rng_x[1] & Time <= rng_x[2] )
                koos <-  koos %>% filter(Time >= rng_x[1] & Time <= rng_x[2] )
            }
      
        } else if (input$range_x != "" &&  input$change_scale == TRUE && vals$Datum) {
            rng_x <- as.Date(strsplit(input$range_x,",")[[1]])
            # rng_x <- c(as.Date(rng_x[1]),as.Date(rng_x[2]))  
        }
      
        #Autoscale if rangeis NOT specified
        else if (input$range_x == "" ||  input$change_scale == FALSE) {
            rng_x <- c(min(klaas$Time), max(klaas$Time))
            #     observe({ print(rng_x) })
        }
    
    
        #Increase rng_x[2] if labels are added BUT NOT for small multiple
        if (input$show_labels_y == TRUE && input$multiples == FALSE) {
            rng_x[2] <- (rng_x[2]-rng_x[1])*.15+rng_x[2]
        }
    }
    else {
        if (input$range_x_2nd != "" &&  input$add_2nd_scale == TRUE && !vals$Datum) {
            rng_x <- as.numeric(strsplit(input$range_x_2nd,",")[[1]])
            observe({ print(rng_x) })
      
      
            #If min>max invert the axis
            if (rng_x[1]>rng_x[2]) {
                p <- p+ scale_x_reverse()
                klaas <-  klaas %>% filter(Time >= rng_x[2] & Time <= rng_x[1] )
                koos <-  koos %>% filter(Time >= rng_x[2] & Time <= rng_x[1] )
            } else {
            #Select timepoints within the set limits
                klaas <-  klaas %>% filter(Time >= rng_x[1] & Time <= rng_x[2] )
                koos <-  koos %>% filter(Time >= rng_x[1] & Time <= rng_x[2] )
            }
      
        } else if (input$range_x_2nd != "" &&  input$add_2nd_scale == TRUE && vals$Datum) {
            rng_x <- as.Date(strsplit(input$range_x_2nd,",")[[1]])
            # rng_x <- c(as.Date(rng_x[1]),as.Date(rng_x[2]))  
        }
      
        #Autoscale if rangeis NOT specified
        else if (input$range_x_2nd == "" ||  input$add_2nd_scale == FALSE) {
            rng_x <- c(min(klaas$Time), max(klaas$Time))
            #     observe({ print(rng_x) })
        }
    
    
        #Increase rng_x[2] if labels are added BUT NOT for small multiple
        if (input$show_labels_y == TRUE && input$multiples == FALSE) {
            rng_x[2] <- (rng_x[2]-rng_x[1])*.15+rng_x[2]
        }
 
    }
    
    
    #Define how colors are used
    klaas <- klaas %>% mutate(id = as.factor(id), unique_id = as.character(unique_id))
    koos <- koos %>% mutate(id = as.factor(id))
    
    number_of_conditions <- nlevels(as.factor(klaas$id))
    if (number_of_conditions == 1) {
      kleur_data <- "unique_id"
    } else if (number_of_conditions > 1) {
        kleur_data <- "id"
      }
  
    if (input$color_data == FALSE) {
      kleur_data <- NULL
    }
    
    if (input$fnt_sz_stim == "") {
      fnt_sz_stim <- 6
    } else {
      fnt_sz_stim <- input$fnt_sz_stim
    }
    
    
    newColors <- NULL
    
    if (input$adjustcolors == 2) {
      newColors <- Tol_bright
    } else if (input$adjustcolors == 3) {
      newColors <- Tol_muted
    } else if (input$adjustcolors == 4) {
      newColors <- Tol_light
    } else if (input$adjustcolors == 6) {
      Okabe_Ito[8] <- line_color
      newColors <- Okabe_Ito
    } else if (input$adjustcolors == 5) {
      newColors <- gsub("\\s","", strsplit(input$user_color_list,",")[[1]])
    }
            
    max_colors <- nlevels(as.factor(klaas$unique_id))
    if(length(newColors) < max_colors) {
      newColors<-rep(newColors,times=(round(max_colors/length(newColors)))+1)
    }
           
    ########## Define how color is mapped onto the data
    #    observe({ print(class(input$colour_list)) })
    if (input$color_stats == FALSE) {
      kleur_stats <- NULL
    } else if (input$color_stats == TRUE) {
      kleur_stats <- "id"
    } 
    
    #### Command to prepare the plot ####
    p <- ggplot(data=klaas, aes_string(x="Time")) 
    
         
    #### plot individual measurements ####
    
    if (input$thicken =="TRUE") {
      multiplier <- 4
    } else if (input$thicken =="FALSE"){
      multiplier <- 1
    }

        
    if (input$data_form == "dataasline") {
      p <- p+ geom_line(data=klaas, aes_string(x="Time", y="Value", group="unique_id", color=kleur_data), size=0.5*multiplier, alpha=input$alphaInput)
      if (!input$color_data) {p <- p+geom_line(aes_string(x="Time", y="Value", group="unique_id"), color=line_color, size=0.5*multiplier, alpha=input$alphaInput)}
    } else if (input$data_form == "dataasdot") {
      p <- p + geom_point(data=klaas, aes_string(x="Time", y="Value", group="unique_id", color=kleur_data), size=1*multiplier, alpha=input$alphaInput)
      if (!input$color_data) {p <- p+geom_point(aes_string(x="Time", y="Value", group="unique_id"), color=line_color, size=1*multiplier, alpha=input$alphaInput)}
    } 

    #### plot stats ####
    
    # if (input$summaryInput == TRUE  && input$add_CI == FALSE) {
    #   p <- p + geom_line(data=koos, aes_string(x="Time", y="mean", group="id", color=kleur_stats),size=2,alpha=input$alphaInput_summ)
    #   if (!input$color_stats) {p <- p + geom_line(data=koos, aes_string(x="Time", y="mean", group="id"),color=line_color,size=2,alpha=input$alphaInput_summ)}
    # } else if (input$summaryInput == TRUE  && input$add_CI == TRUE) {
    #   p <- p + geom_ribbon(data=koos, aes_string(x="Time", ymin="ci_lo", ymax="ci_hi", group="id", fill=kleur_stats), alpha=input$alphaInput_summ/2)
    #   p <- p + geom_line(data=koos, aes_string(x="Time", y="mean", group="id", color=kleur_stats),size=2,alpha=input$alphaInput_summ)
    #   p <- p + guides(fill = FALSE)
    # } else if (input$summaryInput == FALSE  && input$add_CI == TRUE) {
    #   p <- p + geom_ribbon(data=koos, aes_string(x="Time", ymin="ci_lo", ymax="ci_hi", group="id", fill=kleur_stats), alpha=input$alphaInput_summ/2)
    #   if (!input$color_stats) {p <- p + geom_ribbon(data=koos, aes_string(x="Time", ymin="ci_lo", ymax="ci_hi", group="id"), fill=line_color, alpha=input$alphaInput_summ/2)}
    #   
    # }
    
    
    if (input$summaryInput == TRUE) {
      p <- p + geom_line(data=koos, aes_string(x="Time", y="mean", group="id", color=kleur_stats),size=2,alpha=input$alphaInput_summ)
      if (!input$color_stats) {p <- p + geom_line(data=koos, aes_string(x="Time", y="mean", group="id"),color=line_color,size=2,alpha=input$alphaInput_summ)}
    }
   if (input$add_CI == TRUE) {
      p <- p + geom_ribbon(data=koos, aes_string(x="Time", ymin="ci_lo", ymax="ci_hi", group="id", fill=kleur_stats), alpha=input$alphaInput_summ/2)
      if (!input$color_stats) {p <- p + geom_ribbon(data=koos, aes_string(x="Time", ymin="ci_lo", ymax="ci_hi", group="id"), fill=line_color, alpha=input$alphaInput_summ/2)}
      
    }
    
    
    
    
    
    # This needs to go here (before annotations)
    p <- p+ theme_light(base_size = 16)
    if (input$dark) {p <- p+ theme_light_dark_bg(base_size = 16)}
    

    ############## Adjust Y-scaling if necessary ##########
    
    
    #Adjust scale if range for y (min,max) is specified
    if (input$range_y != "" &&  input$change_scale == TRUE) {
      rng_y <- as.numeric(strsplit(input$range_y,",")[[1]])

      #If min>max invert the axis
      if (rng_y[1]>rng_y[2]) {p <- p+ scale_y_reverse()}
      
      upper_y <- rng_y[2]
      lower_y <- rng_y[1]
      

      #Autoscale if rangeis NOT specified
    } else if (input$range_y == "" ||  input$change_scale == FALSE) {

      
      #Read out the current range, this is necessary for annotation of stimulus
      rng_y <- c(NULL,NULL)
      upper_y <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
      lower_y <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]
      
    }
    range_y <- upper_y-lower_y

    #################### Add labels for perturbations #####################
    
    rang <- as.numeric(strsplit(input$stim_range,",")[[1]])
    
    stimText <- c("","","","","")
    
    if (input$indicate_stim == TRUE && input$stim_text !="") {
      stimText <- strsplit(input$stim_text,",")[[1]]
    }

    if (input$indicate_stim == TRUE && input$stim_colors !="") {
      stimColors <- gsub("\\s","", strsplit(input$stim_colors,",")[[1]])
      
    } else if (input$indicate_stim == TRUE && input$stim_colors =="") {
      stimColors <- "black"
    }
    
    # if a stimulus is applied
    if (input$indicate_stim == TRUE) {
      
      p <- p  +  theme(plot.margin = unit(c(3,1,1,1), "lines"))
      p <- p + coord_cartesian(xlim=c(rng_x[1],rng_x[2]),ylim=c(lower_y,upper_y),clip = 'off')
      
      #If only one number is entered, a vertical line is added
      if (length(rang) ==1) {
        p <- p + geom_vline(xintercept=rang[1], black="orange", size=1)
      }
      
      nsteps = floor(length(rang)/2)
      #Repeat the colors if needed
      if(length(stimColors) < nsteps) {
        stimColors<-rep(stimColors,times=(round(nsteps/length(stimColors)))+1)
      }
      
      
      if(input$stim_shape == "bar") {
          for (i in 0:nsteps) {
             p <- p + annotate("rect", xmin=rang[(i*2)+1], xmax=rang[(i*2)+2], ymin=upper_y+.02*range_y, ymax=upper_y+.05*range_y, alpha=0.8, fill=stimColors[i+1])
             p <- p + annotate("text", x=rang[(i*2)+1]+0.5*(rang[(i*2)+2]-rang[(i*2)+1]), y=upper_y+.1*range_y, alpha=1, color=stimColors[i+1], size=fnt_sz_stim,label=paste(stimText[i+1]))
          }
      } else if (input$stim_shape == "box") {
        
          for (i in 0:nsteps) {
              p <- p + annotate("rect", xmin=rang[(i*2)+1], xmax=rang[(i*2)+2], ymin=-Inf, ymax=Inf, alpha=0.1, fill=stimColors[i+1])
              p <- p + annotate("text", x=rang[(i*2)+1]+0.5*(rang[(i*2)+2]-rang[(i*2)+1]), y=upper_y+.1*range_y, alpha=1, color=stimColors[i+1], size=fnt_sz_stim,label=paste(stimText[i+1]))
          }
          
      } else if (input$stim_shape == "both") {
          for (i in 0:nsteps) {
            p <- p + annotate("rect", xmin=rang[(i*2)+1], xmax=rang[(i*2)+2], ymin=-Inf, ymax=Inf, alpha=0.1, fill=stimColors[i+1])
            p <- p + annotate("rect", xmin=rang[(i*2)+1], xmax=rang[(i*2)+2], ymin=upper_y+.02*range_y, ymax=upper_y+.05*range_y, alpha=0.8, fill=stimColors[i+1])
            p <- p + annotate("text", x=rang[(i*2)+1]+0.5*(rang[(i*2)+2]-rang[(i*2)+1]), y=upper_y+.1*range_y, alpha=1, color=stimColors[i+1], size=fnt_sz_stim,label=paste(stimText[i+1]))
          }
        
      }
    } else {
      p <- p + coord_cartesian(xlim=c(rng_x[1],rng_x[2]),ylim=c(lower_y,upper_y))
    }
    

    
########## Do some formatting of the lay-out ##########
    


    # if title specified
    if (input$add_title == TRUE) {
      #Add line break to generate some space
      title <- paste(input$title, "\n",sep="")
      p <- p + labs(title = title)
    }

    # # if labels specified
    if (input$label_axes)
      p <- p + labs(x = input$lab_x, y = input$lab_y)

    

    # # if font size is adjusted
    if (input$adj_fnt_sz) {
      p <- p + theme(axis.text = element_text(size=input$fnt_sz_ax))
      p <- p + theme(axis.title = element_text(size=input$fnt_sz_labs))
      p <- p + theme(plot.title = element_text(size=input$fnt_sz_title))
      }
    
    #remove legend (if selected)
    if (input$add_legend == FALSE) {  
      p <- p + theme(legend.position="none")
    }
    

    ################################ Add labels  ####################
    # For traces in case of one conditions
    # For stats in case of multiple
     
    
    
    
    #Generate a dataframe with the labels
    if (input$show_labels_y == TRUE && input$multiples == FALSE) {
            
            # If summary is selected, label mean trace
            if (input$summaryInput == TRUE) {
              df_label <- koos %>% group_by(id) %>% filter(Time==last(Time))
              
            #If  summary is not selected, label all traces  
            } else if (input$summaryInput == FALSE) {
              df_label <- klaas %>% filter(Time==last(Time))
              df_label$unique_id <- gsub(".*_","",df_label$unique_id)
            }
      
 
            
            if (input$summaryInput == FALSE && input$color_data == FALSE) {
              p <- p + geom_label_repel(data = df_label, aes_string(label='unique_id', x='Time', y='Value'),
                                        fill = 'black',
                                        fontface = 'bold', color = 'white', size=6,
                                        nudge_x      = 20,
                                        direction    = "y",
                                        hjust        = 0,
                                        point.padding = unit(1, 'lines'),
                                        segment.color = 'grey50',
                                        #xlim = c(input$lim_start_x_1, input$lim_end_x_1),
                                        segment.size = 0.5)
              

            } else  if (input$summaryInput== FALSE && input$color_data == TRUE) {
                #  for multiple conditions fill = id, but for a single condition use fill='unique_id'
                p <- p + geom_label_repel(data = df_label, aes_string(label='unique_id', x='Time', y='Value', fill=kleur_data),

                                          
                                          
                                          fontface = 'bold', color = 'white', size=6,
                                          nudge_x      = 20,
                                          direction    = "y",
                                          hjust        = 0,
                                          point.padding = unit(1, 'lines'),
                                          segment.color = 'grey50',
                                          #  xlim = c(input$lim_start_x_1, input$lim_end_x_1),
                                          segment.size = 0.5)
                                          
                
            } else  if (input$summaryInput== TRUE && input$color_stats == FALSE) {

              p <- p + geom_label_repel(data = df_label, aes_string(label='id', x='Time', y='mean'),
                                        fill = 'black',
                                        fontface = 'bold', color = 'white', size=6,
                                        nudge_x      = 20,
                                        direction    = "y",
                                        hjust        = 0,
                                        point.padding = unit(1, 'lines'),
                                        segmaent.color = 'grey50',
                                        #xlim = c(input$lim_start_x_1, input$lim_end_x_1),
                                        segment.size = 0.5)
              
              
              
            } else  if (input$summaryInput== TRUE && input$color_stats == TRUE) {
              
              p <- p + geom_label_repel(data = df_label, aes_string(label='id', x='Time', y='mean', fill='id'),
                                        
                                        fontface = 'bold', color = 'white', size=6,
                                        nudge_x      = 20,
                                        direction    = "y",
                                        hjust        = 0,
                                        point.padding = unit(1, 'lines'),
                                        segment.color = 'grey50',
                                        #xlim = c(input$lim_start_x_1, input$lim_end_x_1),
                                        segment.size = 0.5)
              
            }
      }


    
    if (input$show_labels_y == TRUE && input$multiples == TRUE) {
    #Show labels in upper right corner
      
      
      if(number_of_conditions == 1) {

        #show unique_id in upper right corner
        df_label <- klaas %>% filter(Time==last(Time)) 
        
      } else if (number_of_conditions > 1) {

        #show id in upper right corner
        df_label <- koos %>% group_by(id) %>% filter(Time==last(Time)) %>% mutate(unique_id=id)
      }
      
      if (number_of_conditions == 1 && input$color_data == FALSE) {

        # p <- p + geom_label(data = df_label, aes(label=unique_id,x=Inf,y=Inf),
        # I replaced the line above for better support of object labeling (avoiding long, ugly names)
        
        p <- p + geom_label(data = df_label, aes(label=Sample,x=Inf,y=Inf),
                                  fill = 'black',
                                  fontface = 'bold', color = 'white', size=5,
                                  vjust = 1,
                                  hjust = 1)
       
        
      } else if (number_of_conditions == 1 && input$color_data == TRUE) {
        p <- p + geom_label(data = df_label, aes(label=Sample,x=Inf,y=Inf, fill=Sample),
                            fontface = 'bold', color = 'white', size=5,
                            vjust = 1,
                            hjust = 1)
        
      } else  if (number_of_conditions > 1 && input$color_data == FALSE) {
        
        p <- p + geom_label_repel(data = df_label, aes_string(label='id', x=Inf, y=Inf),
                             fill = 'black',
                             fontface = 'bold', color = 'white', size=8,
                             vjust = 1,
                             hjust = 1)
        
      } else  if (number_of_conditions > 1 && input$color_data == TRUE) {
        
        p <- p + geom_label_repel(data = df_label, aes_string(label='id', x=Inf, y=Inf, fill='id'),
                             fontface = 'bold', color = 'white', size=8,
                             vjust = 1,
                             hjust = 1)
        
      }
    }
    
    #Add legend if specified
    if (input$add_legend == TRUE) {
       p <- p + labs(color = input$legend_title, fill=input$legend_title)
    }
    
    # TODO: check!!
    if (input$axes_limits)
      p <- p + coord_cartesian(xlim = c(input$lim_start_x_1, input$lim_end_x_1))
      #p <- p + xlim(input$lim_start_x_1, input$lim_end_x_1)
      #p <- p + scale_x_continuous(limits = c(input$lim_start_x_1, input$lim_end_x_1))



    # if log-scale checked specified
    if (input$scale_log_10)
      p <- p + scale_y_log10() 
    
    
    #remove gridlines (if selected)
    if (input$no_grid == TRUE) {  
      p <- p+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())
    }

    if (input$adjustcolors >1 && input$adjustcolors < 7) {
     p <- p+ scale_color_manual(values=newColors)
     p <- p+ scale_fill_manual(values=newColors)
    } else if (input$adjustcolors ==7) {
      p <- p+ scale_colour_viridis_d()
      p <- p+ scale_fill_viridis_d()      
    }
    
    if (input$multiples == TRUE) {
      if (number_of_conditions == 1) {
                  p <- p+ facet_wrap(~unique_id)
                  #Remove the strip above the individual panels
                  p <- p + theme(strip.background = element_blank(), strip.text = element_blank(), panel.spacing.y = unit(.9, "lines"),panel.spacing.x = unit(.6, "lines"))
      } else if (number_of_conditions > 1) {
                    p <- p+ facet_grid(id~.)
      }
                    
    }
    
    # Hide the facet 'strips' when object are labeled
    if (input$show_labels_y == TRUE && number_of_conditions > 1) {
      p <- p + theme(strip.background = element_blank(), strip.text = element_blank(), panel.spacing.y = unit(.5, "lines"),panel.spacing.x = unit(.5, "lines"))
      
    }

    
        #Extend the canvas to avoid clipping of x-axis labels 
        p <- p + theme(plot.margin = unit(c(5.5,16,5.5,5.5), "pt"))
        
        #Remove upper and right axis
        p <- p + theme(panel.border = element_blank())
        p <- p + theme(axis.line.x  = element_line(colour = line_color), axis.line.y  = element_line(colour = line_color))
        

        
    return(p)
    
    return(p)
}
