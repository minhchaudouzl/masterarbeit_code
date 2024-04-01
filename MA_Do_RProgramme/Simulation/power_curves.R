library(rlist)
library(ggplot2)

mypath <- "D:/MA_Do_RProgramme"

# Power- und Effektmaßdaten zusammenführen
format_power_results <- function(power_results){
  all_results <- as.data.frame(do.call(rbind, power_results))
  all_results <- data.frame(power = unlist(all_results[, 1]),
                            effect = unlist(all_results[, 2]))
  
  return(all_results)
}

gleason_log <- list.load(paste0(mypath, "/Simulation/Gleason/gleason_log.rds"))
gleason_rmst <- list.load(paste0(mypath, "/Simulation/Gleason/gleason_rmst.rds"))
gleason_yp <- list.load(paste0(mypath, "/Simulation/Gleason/gleason_yp.rds"))


jorgensen_log <- list.load(paste0(mypath, "/Simulation/Jorgensen/jorgensen_log.rds"))
jorgensen_rmst <- list.load(paste0(mypath, "/Simulation/Jorgensen/jorgensen_rmst.rds"))
jorgensen_yp <- list.load(paste0(mypath, "/Simulation/Jorgensen/jorgensen_yp.rds"))


leon_log <- list.load(paste0(mypath, "/Simulation/Leon/leon_log.rds"))
leon_rmst <- list.load(paste0(mypath, "/Simulation/Leon/leon_rmst.rds"))

leon_yp <- list.load(paste0(mypath, "/Simulation/Leon/leon_yp.rds"))
# Short-term HR größer als 1.5 entfernen
rows <- row.names(format_power_results(leon_yp))
leon_yp_indices <- as.numeric(rows[which(format_power_results(leon_yp)$effect>1.5)])
leon_yp <- leon_yp[-leon_yp_indices]


makkar2012_log <- list.load(paste0(mypath, "/Simulation/Makkar2012/makkar2012_log.rds"))
makkar2012_rmst <- list.load(paste0(mypath, "/Simulation/Makkar2012/makkar2012_rmst.rds"))

makkar2012_yp <- list.load(paste0(mypath, "/Simulation/Makkar2012/makkar2012_yp.rds"))
# Short-term HR größer als 10 entfernen
rows <- row.names(format_power_results(makkar2012_yp))
makkar2012_yp_indices <- as.numeric(rows[which(format_power_results(makkar2012_yp)$effect>10)])
makkar2012_yp <- makkar2012_yp[-makkar2012_yp_indices]


makkar2020_log <- list.load(paste0(mypath, "/Simulation/Makkar2020/makkar2020_log.rds"))
makkar2020_rmst <- list.load(paste0(mypath, "/Simulation/Makkar2020/makkar2020_rmst.rds"))

makkar2020_yp <- list.load(paste0(mypath, "/Simulation/Makkar2020/makkar2020_yp.rds"))
# Short-term HR größer als 1.3 entfernen
rows <- row.names(format_power_results(makkar2020_yp))
makkar2020_yp_indices <- as.numeric(rows[which(format_power_results(makkar2020_yp)$effect>1.3)])
makkar2020_yp <- makkar2020_yp[-makkar2020_yp_indices]


mack_log <- list.load(paste0(mypath, "/Simulation/Mack/mack_log.rds"))
mack_rmst <- list.load(paste0(mypath, "/Simulation/Mack/mack_rmst.rds"))

mack_yp <- list.load(paste0(mypath, "/Simulation/Mack/mack_yp.rds"))
rows <- row.names(format_power_results(mack_yp))
mack_yp_indices <- as.numeric(rows[which(format_power_results(mack_yp)$effect>1.25)])
mack_yp <- mack_yp[-mack_yp_indices]


popma_log <- list.load(paste0(mypath, "/Simulation/Popma/popma_log.rds"))
popma_rmst <- list.load(paste0(mypath, "/Simulation/Popma/popma_rmst.rds"))
popma_yp <- list.load(paste0(mypath, "/Simulation/Popma/popma_yp.rds"))


thyregod_log <- list.load(paste0(mypath, "/Simulation/Thyregod/thyregod_log.rds"))

thyregod_rmst <- list.load(paste0(mypath, "/Simulation/Thyregod/thyregod_rmst.rds"))
# Ein Ausreißer wurde entfernt
rows <- row.names(format_power_results(thyregod_rmst))
thyregod_rmst_indices <- as.numeric(rows[which(format_power_results(thyregod_rmst)$effect< -2)])
thyregod_rmst <- thyregod_rmst[-thyregod_rmst_indices]

thyregod_yp <- list.load(paste0(mypath, "/Simulation/Thyregod/thyregod_yp.rds"))


all_log <- list.load(paste0(mypath, "/Simulation/all/all_log.rds"))
rows <- row.names(format_power_results(all_log))
all_log_indices <- as.numeric(rows[which(format_power_results(all_log)$effect>2.5)])
all_log <- all_log[-all_log_indices]

all_rmst <- list.load(paste0(mypath, "/Simulation/all/all_rmst.rds"))
rows <- row.names(format_power_results(all_rmst))
all_rmst_indices <- as.numeric(rows[which(format_power_results(all_rmst)$effect>20 |
                                          format_power_results(all_rmst)$effect< -30)])
all_rmst <- all_rmst[-all_rmst_indices]

all_yp <- list.load(paste0(mypath, "/Simulation/all/all_yp.rds"))
rows <- row.names(format_power_results(all_yp))
all_yp_indices <- as.numeric(rows[which(format_power_results(all_yp)$effect>2)])
all_yp <- all_yp[-all_yp_indices]

all_indices <- union(union(all_log_indices, all_rmst_indices), all_yp_indices)



output_dir <- paste0(mypath, "/Simulation")

plot_power_curve_smooth <- function(power_results_log, power_results_rmst, 
                                    power_results_yp, study) {
  # Plotten der Powerkurve
  log_results <- format_power_results(power_results_log)
  
  ggplot(log_results, aes(x = effect, y = power)) +
    geom_smooth(method = 'loess', se = TRUE, color = 'blue') +
    geom_point() + 
    labs(x = 'Hazard Ratio', y = 'Power') +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    scale_y_continuous(limits = c(0, 1)) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm")) +
    geom_hline(yintercept = 0.05, linetype = "dashed")
  
  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_log_smooth_", study, ".png")), 
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)


  rmst_results <- format_power_results(power_results_rmst)
  ggplot(rmst_results, aes(x = effect, y = power)) +
    geom_smooth(method = 'loess',se = TRUE, color = 'red') +
    geom_point() +
    labs(x = 'RMST-Differenz', y = 'Power') +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    scale_y_continuous(limits = c(0, 1)) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm")) +
    geom_hline(yintercept = 0.05, linetype = "dashed")

  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_rmst_smooth_", study, ".png")),
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)


  yp_results <- format_power_results(power_results_yp)
  ggplot(yp_results, aes(x = effect, y = power)) +
    geom_smooth(method = 'loess', se = TRUE, color = 'purple') +
    geom_point() +
    labs(x = 'Short-term Hazard Ratio', y = 'Power') +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    scale_y_continuous(limits = c(0, 1)) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm")) +
    geom_hline(yintercept = 0.05, linetype = "dashed")

  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_yp_smooth_", study, ".png")),
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)
}

# Powerkurve plotten
plot_power_curve_smooth(gleason_log, gleason_rmst, gleason_yp, "Gleason")
plot_power_curve_smooth(jorgensen_log, jorgensen_rmst, jorgensen_yp, "Jorgensen")
plot_power_curve_smooth(leon_log, leon_rmst, leon_yp, "Leon")
plot_power_curve_smooth(mack_log, mack_rmst, mack_yp, "Mack")
plot_power_curve_smooth(makkar2012_log, makkar2012_rmst, makkar2012_yp, "Makkar2012")
plot_power_curve_smooth(makkar2020_log, makkar2020_rmst, makkar2020_yp, "Makkar2020")
plot_power_curve_smooth(popma_log, popma_rmst, popma_yp, "Popma")
plot_power_curve_smooth(thyregod_log, thyregod_rmst, thyregod_yp, "Thyregod")
plot_power_curve_smooth(all_log, all_rmst, all_yp, "all")


# Powern anderer Effektmaße über einem jeden Effektmaß zeichnen
power_plot_combined <- function(power_results_log, power_results_rmst, 
                                power_results_yp, study){
  # Plotten der Powerkurve
  log_results <- format_power_results(power_results_log)
  colnames(log_results) <- c("power_lr", "hr")
  rmst_results <- format_power_results(power_results_rmst)
  colnames(rmst_results) <- c("power_rmst", "rmst")
  yp_results <- format_power_results(power_results_yp)
  colnames(yp_results) <- c("power_yp", "short_term_HR")
  
  results <- cbind(cbind(log_results, rmst_results), yp_results)

  
  ggplot() +
    geom_point(data = results, aes(x = hr, y = power_lr), color = "dodgerblue") +
    geom_point(data = results, aes(x = hr, y = power_rmst), color = "indianred1") +
    geom_point(data = results, aes(x = hr, y = power_yp), color = "darkorchid3") +
    geom_smooth(data = results, aes(x = hr, y = power_lr), 
                method = 'loess', se = FALSE, color = "blue") +
    geom_smooth(data = results, aes(x = hr, y = power_rmst), 
                method = 'loess', se = FALSE, color = "red") +
    geom_smooth(data = results, aes(x = hr, y = power_yp), 
                method = 'loess', se = FALSE, color = "darkviolet") +
    labs(x = "Hazard Ratio", y = "Power") +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    ylim(0, 1) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm"))
  
  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_log_combined_", study, ".png")), 
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)
  
  ggplot() +
    geom_point(data = results, aes(x = rmst, y = power_lr), color = "dodgerblue") +
    geom_point(data = results, aes(x = rmst, y = power_rmst), color = "indianred1") +
    geom_point(data = results, aes(x = rmst, y = power_yp), color = "darkorchid3") +
    geom_smooth(data = results, aes(x = rmst, y = power_lr), 
                method = 'loess', se = FALSE, color = "blue") +
    geom_smooth(data = results, aes(x = rmst, y = power_rmst), 
                method = 'loess', se = FALSE, color = "red") +
    geom_smooth(data = results, aes(x = rmst, y = power_yp), 
                method = 'loess', se = FALSE, color = "darkviolet") +
    labs(x = "RMST-Differenz", y = "Power") +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    ylim(0, 1) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm"))
  
  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_rmst_combined_", study, ".png")), 
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)
  
  ggplot() +
    geom_point(data = results, aes(x = short_term_HR, y = power_lr), color = "dodgerblue") +
    geom_point(data = results, aes(x = short_term_HR, y = power_rmst), color = "indianred1") +
    geom_point(data = results, aes(x = short_term_HR, y = power_yp), color = "darkorchid3") +
    geom_smooth(data = results, aes(x = short_term_HR, y = power_lr), 
                method = 'loess', se = FALSE, color = "blue") +
    geom_smooth(data = results, aes(x = short_term_HR, y = power_rmst), 
                method = 'loess', se = FALSE, color = "red") +
    geom_smooth(data = results, aes(x = short_term_HR, y = power_yp), 
                method = 'loess', se = FALSE, color = "darkviolet") +
    labs(x = "Short-term Hazard Ratio", y = "Power") +
    theme(axis.text = element_text(size = rel(1.8)), axis.title = element_text(size = rel(1.8))) +
    ylim(0, 1) + theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm"))
  
  ggsave(file.path(paste0(output_dir, "/", study), paste0("power_yp_combined_", study, ".png")), 
         units = "px", height=480*(300/72), width=480*(300/72), dpi = 300)
}


power_plot_combined(gleason_log, gleason_rmst, gleason_yp, "Gleason")  
power_plot_combined(jorgensen_log, jorgensen_rmst, jorgensen_yp, "Jorgensen")
power_plot_combined(leon_log[-leon_yp_indices], leon_rmst[-leon_yp_indices], leon_yp, "Leon")
power_plot_combined(mack_log[-mack_yp_indices], mack_rmst[-mack_yp_indices], mack_yp, "Mack")
power_plot_combined(makkar2012_log[-makkar2012_yp_indices], 
                    makkar2012_rmst[-makkar2012_yp_indices], makkar2012_yp, "Makkar2012")
power_plot_combined(makkar2020_log[-makkar2020_yp_indices], 
                    makkar2020_rmst[-makkar2020_yp_indices], makkar2020_yp, "Makkar2020")
power_plot_combined(popma_log, popma_rmst, popma_yp, "Popma")
power_plot_combined(thyregod_log[-thyregod_rmst_indices], thyregod_rmst, 
                    thyregod_yp[-thyregod_rmst_indices], "Thyregod")

# Indices von allen drei Tests rausnehmen
all_log <- list.load(paste0(mypath, "/Simulation/all/all_log.rds"))
all_rmst <- list.load(paste0(mypath, "/Simulation/all/all_rmst.rds"))
all_yp <- list.load(paste0(mypath, "/Simulation/all/all_yp.rds"))
power_plot_combined(all_log[-all_indices], all_rmst[-all_indices], all_yp[-all_indices], "all")
  