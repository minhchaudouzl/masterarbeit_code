library(survivalnma)

# Hier muss der entsprechende Pfad zum Ordner mit den Rekonstruktionsdaten eingefügt werden
setwd("D:/MA_Do_RProgramme/Rekonstruktion/Jorgensen")

data_blau <- read.csv2("blau.csv", sep = ";", header = FALSE)
data_rot <- read.csv2("rot.csv", sep = ";", header = FALSE)

# Spalten sinnvoll umbenennen
colnames(data_blau) <- c("T_k", "S_k")
colnames(data_rot) <- c("T_k", "S_k")

# T_k in Monate umrechnen
data_blau$T_k <- data_blau$T_k * 12
data_rot$T_k <- data_rot$T_k * 12


# scale your y to be in [0, 1] range
# Werte außerhalb des Bereichs [0, 100] raus
data_blau <- subset(data_blau, S_k >= 0 & S_k <= 100)
data_rot <- subset(data_rot, S_k >= 0 & S_k <= 100)

# guyot.method(data_blau$T_k, data_blau$S_k, nrisk_blau$trisk, nrisk_blau$nrisk)


# mehrfach vorkommende T_k Werte rausnehmen
data_blau <- data_blau[!duplicated(data_blau$T_k), ]
data_rot <- data_rot[!duplicated(data_rot$T_k), ]





#' sort_T_k gibt einen veraederten Dataframe zurueck, indem alle T_k Werte, die 
#' nicht der Groeße nach sortiert sind, entfernt werden
#' 
#' @param df    Ein nach der y-Achse sortierter Dataframe
#'
#' @return df_filtered   Dataframe 


sort_T_k <- function(df){
  # Erstelle eine leere Liste, um die Zeilen zu speichern
  filtered_rows <- list()
  
  # Initialisiere den letzten behaltenen Wert
  last_kept_value <- -Inf
  
  # Durchlaufe die Zeilen des Data Frames
  for (i in 1:nrow(df)) {
    current_value <- df$T_k[i]
    if (current_value >= last_kept_value) {
      filtered_rows[[length(filtered_rows) + 1]] <- df[i, ]
      last_kept_value <- current_value
    }
  }
  
  # Erstelle einen Data Frame aus den gespeicherten Zeilen
  df_filtered <- do.call(rbind, filtered_rows)
  
  # Zurückgebliebenes DataFrame anzeigen
  return(df_filtered)
}


entferne <- function(df, x){
  ent <- c()
  for (i in 1:(length(x)-1)) {
    if(diff(x,1)[i] < 0){
      ent <- c(ent, (i+1))
    }
  }
  ent
  df <- df[-ent,]
  return(df)
}



### Sortierung der roten Daten
# TAVR
data_rot <- data.frame(c(0, data_rot$T_k), c(100, data_rot$S_k))

colnames(data_rot) <- c("T_k", "S_k")
# nach S_k
data_rot <- data_rot[order(data_rot$S_k, decreasing = TRUE), ]
# nach T_k
data_rot <- sort_T_k(data_rot)

# Überprüfung
data_rot$T_k == sort(data_rot$T_k)

data_rot$T_k <- as.numeric(data_rot$T_k)
data_rot$S_k <- as.numeric(data_rot$S_k)

# T_k auf range(0,1) bringen
data_rot$S_k <- data_rot$S_k / 100
write.csv2(data_rot, file = "jorgensen_tavr.csv", row.names = FALSE)

tot.events_data_rot <- "NA"
arm.id_data_rot <- 1


### Sortierung der blauen Daten
# SAVR
data_blau <- data.frame(c(0, data_blau$T_k), c(100, data_blau$S_k))

colnames(data_blau) <- c("T_k", "S_k")
# nach S_k
data_blau <- data_blau[order(data_blau$S_k, decreasing = TRUE), ]
# nach T_k
data_blau <- sort_T_k(data_blau)

# Überprüfung
data_blau$T_k == sort(data_blau$T_k)

data_blau$T_k <- as.numeric(data_blau$T_k)
data_blau$S_k <- as.numeric(data_blau$S_k)
# T_k auf range(0,1) bringen
data_blau$S_k <- data_blau$S_k / 100
write.csv2(data_blau, file = "jorgensen_savr.csv", row.names = FALSE)

tot.events_data_blau <- "NA"
arm.id_data_blau <- 0



nrisk_blau <- read.csv2("nrisk_blau.csv", sep = ";")
nrisk_rot <- read.csv2("nrisk_rot.csv", sep = ";")







#' guyot rekonstruiert Daten aus Daten, die aus z.B. Kaplan-Meier-Kurven 
#' abgelesen wurden 
#' 
#' @param dataframe     Ein zweispaltiger Dataframe, der nach der y-Achse 
#'                      sortierte Wertepaare enthält, wobei alle x-Werte auch 
#'                      aufsteigend sortiert sind
#' @param risktable     Eine Tabelle (nrisk), wie im Paper von guyot beschrieben 
#' @param tot.events    Gesamtanzahl der Ereignisse der Gruppe (auch NA moeglich)
#' @param armID         zu welcher Gruppe gehoeren die Daten (wichtig, wenn man 
#'                      2 verschiedene Gruppen hinterher vergleichen moechte)
#'
#' @return              Matrix mit Eintraegen, wann ein Event stattgefunden hat

guyot <- function(dataframe, risktable, tot.events, armID, IPD){
  t.S <- dataframe[,1]
  S <- dataframe[,2]
  
  t.risk <- risktable[,1]
  lower <- risktable[,3]
  upper <- risktable[,4]
  n.risk <- risktable[,2]
  n.int <- length(n.risk)
  n.t <- upper[n.int]
  
  arm <- rep(armID, n.risk[1])
  n.censor <- rep(0, (n.int - 1))
  n.hat <- rep(n.risk[1] + 1, n.t)
  cen <- rep(0, n.t)
  d <-  rep(0, n.t)
  KM.hat <- rep(1, n.t)
  last.i <- rep(1, n.int)
  sumdL <- 0
  
  # Annahme censoring ist unbekannt, aber constant
  # erste Vermutung 
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
      # Erste censoring Annahme, wenn keine zensierten Daten da sind, dann werden 
      # die unter Risiko stehenden Personen und die Wahrscheinlichkeiten am 
      # Anfang/Ende des Intervalls zu überleben miteinander verrechnet
      n.censor[i] <- round(n.risk[i]*S[lower[i+1]]/S[lower[i]] - n.risk[i+1])
      # Verteilung der Zensierungen gleichmäßig auf das Intervall 
      while((n.hat[lower[i+1]] > n.risk[i+1])||((n.hat[lower[i+1]] < n.risk[i+1]) && (n.censor[i] > 0))){
        # wenn der Zensierungswert kleiner als 0 ist, bekommt er den Wert 0
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        # alle zensierten Werte größer null werden gleichmäßig verteilt
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            # (4) Formel
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          # Verteilung der Daten gleichmäßig über das Intervall
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],
                                       plot=F)$counts
        }
        # Kalkulation erstellen anhand von Eventanzahl an den KM-Koordinaten und 
        # den unter Risiko stehenden Personen zum nächsten Zeitpunkt
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        # mindestens ein event am Beginn und Ende des Intervalls
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k] <- 0
            KM.hat[k] <- 1
          }
          else {
            # (6) Formel
            d[k] <- round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k] <- KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          # Formel (7)
          n.hat[k+1] <- n.hat[k] - d[k] - cen[k]
          if (d[k] != 0) last <- k
        }
        n.censor[i] <- n.censor[i]+ (n.hat[lower[i+1]] - n.risk[i+1])
      }
      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1] <- n.hat[lower[i+1]]
      last.i[(i+1)] <- last
    }
  }
  
  
  # Zeitintervall n.int.
  if (n.int>1){
    # Annahme: gleiche Zensierungsrate als Durchschnitt über den vorherigen 
    # Intervallen
    n.censor[n.int] <- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-                                                      t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  
  if (n.int==1){
    n.censor[n.int]<-0
  }
  
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  
  if (n.censor[n.int] > 0){
    cen.t <- rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j] <- t.S[lower[n.int]] + j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)] <- hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                               plot=F)$counts
  }
  
  # Wdh. der Schritte
  n.hat[lower[n.int]] <- n.risk[n.int]
  last <- last.i[n.int]
  # Seite 11 im Paper Schritt 2 und 3 nochmal laufen lassen
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  
  # Wenn Anzahl total events gegeben, dann diesser Teil, Schritt 7
  if (tot.events != "NA"){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      # Wenn bereits zu viele events auftreten, dann sollen alle folgenden events 
      # und censoring auf 0 gesetzt werden
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Adjustieren nach Schritt 8
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                                   plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            # Ueberpruefung, dass die unter Risiko stehenden Personen nicht negativ ist
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  
  matrix1 <- matrix(c(t.S,n.hat[1:n.t],d,cen),ncol=4,byrow=F)
  
  # Ab hier erstellen der Tabellen
  
  ### IPD (induvidual patient data) erstellen ###
  # Vektor initialisiereb
  t.IPD<-rep(t.S[n.t],n.risk[1])
  event.IPD<-rep(0,n.risk[1])
  # event time und event indicator (=1) für jedes event schreiben, separierte Spalten in t.IPD und event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  # event time und event indicator (=0) für jede Zensierung schreiben, separierte Spalten in t.IPD und event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }
  # Output IPD
  IPD<-matrix(c(t.IPD,event.IPD,arm),ncol=3,byrow=F)
  
  return(IPD)
  #ifelse(IPD == FALSE,return(matrix1), )
  
}

# number at risk werte manuell?

data_rot_guyot <- guyot(dataframe = data_rot, risktable = nrisk_rot, tot.events = tot.events_data_rot, armID = arm.id_data_rot)

data_blau_guyot <- guyot(dataframe = data_blau, risktable = nrisk_blau, tot.events = tot.events_data_blau, armID = arm.id_data_blau)


write.table(x = data_rot_guyot,
            file = "jorgensen_rot_guyot.csv",
            sep = ";",
            row.names = FALSE,
            col.names = FALSE)

write.table(x = data_blau_guyot,
            file = "jorgensen_blau_guyot.csv",
            sep = ";",
            row.names = FALSE,
            col.names = FALSE)



