# Load necessary packages
packages <- c("data.table", "tidyverse", "keyring", "blastula", "dtplyr", 
              "naniar", "network", "sna", "Matrix", "haven", "xtable", 
              "ndtv", "networkDynamic")
needed_packages <- setdiff(packages, rownames(installed.packages()))
lapply(needed_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

# Upload data from specified paths
e0 <- read.csv("Path to relationship data")
House_16_Presence <- read.csv("path to presence data")

# Assign weights to edges based on friendship ratings
e0$wt <- if_else(e0$frd < 3, 1, 0)
e1 <- e0[, 1:4]

# Helper function to find a specific vertex by name
get_vertex_id_by_name <- function(net, name) {
  return (which(get.vertex.attribute(net, "vertex.names") == name))
}

# Count the number of distinct waves
wave_num <- unique(e0$h)

# Create a vector of wave table names
wave_names <- paste0("w", 1:7) 

# Initialize wave data frames as empty data frames if they don't exist
for (name in wave_names) {
  if (!exists(name)) {
    assign(name, data.frame())
  }
}

# Use mget to retrieve a list of these data frames
wave_tables <- mget(wave_names)

# Create separate tables for each wave
for(wave in wave_num){
  wave_table_name <- paste("w", wave, sep = "")
  assign(wave_table_name, e1 %>% filter(h == wave))
}

# Initialize a network and assign names to each vertex
all_vertices <- unique(c(e0$source, e0$target))
net_dynamic <- network.initialize(length(all_vertices), directed = TRUE)
set.vertex.attribute(net_dynamic, "vertex.names", all_vertices)

# Create tables of valid edges for each wave
for(wave in wave_num){
  valid_table_name <- paste("validEdge_w", wave, sep = "")
  current_table_name <- paste("w", wave, sep = "")
  cur_table <- get(current_table_name)
  valid_edges_table <- cur_table %>% filter(wt == 1)
  assign(valid_table_name, valid_edges_table, envir = .GlobalEnv )
}

# Activates all vertices
activate.vertices(net_dynamic, onset = 0, terminus = length(wave_num) + 1)

# Deactivates the vertices in the first interval
deactivate.vertices(net_dynamic, onset = 0, length = 1)

# Creates separate lists of people not in the house for each wave
for(wave in wave_num){
  cur_wave <- paste("w", wave, sep = "")
  cur_name <- paste("notIn_w", wave, sep = "")
  individual_ids <- House_16_Presence %>%
    filter(!!sym(cur_wave) == 0) %>%
    pull(ID)
  assign(cur_name, individual_ids)
}

# Deactivates vertices for individuals not in the house
for(wave in wave_num){
  in_name <- paste("notIn_w", wave, sep = "")
  cur_active <- get(in_name)
  
  # Loop through the list of people not in the house
  for(i in 1:length(cur_active)){
    cur_id <- get_vertex_id_by_name(net_dynamic, cur_active[[i]])
    deactivate.vertices(net_dynamic, onset = wave, length = 1, v = cur_id, deactivate.edges = FALSE)
  }
}

# Activates valid edges at each wave
for(wave in 1:length(wave_num)){
  #Find the current wave
  valid_edges_table_name <- paste("validEdge_w", wave, sep = "")
  valid_edges <- get(valid_edges_table_name)
  
  
  #Actives edges at current wave
  for(i in 1:nrow(valid_edges)){
    source_vertex <- which(get.vertex.attribute(net_dynamic, "vertex.names") == valid_edges[i, "source"])
    target_vertex <- which(get.vertex.attribute(net_dynamic, "vertex.names") == valid_edges[i, "target"])
    
    #Check if the vertex exist
    if(length(source_vertex) > 0 && length(target_vertex) > 0){
      
      #Check if the edge already existed, adds new edge and activates if not
      existing_edge_ids <- get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex)
      if(length(existing_edge_ids) == 0){
        add.edge(net_dynamic, tail = source_vertex, head = target_vertex, activate = FALSE)
        activate.edges(net_dynamic, onset = wave, terminus = wave,
                       e = get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex))
      }
      
      #Activates an existing edge
      activate.edges(net_dynamic, onset = wave, terminus = wave,
                     e = get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex))
    }
  }
}

# Initializes network with vertices
all_vertices <- unique(c(e0$source, e0$target))
net_dynamic <- network.initialize(length(all_vertices), directed = TRUE)
network::set.vertex.attribute(net_dynamic, "vertex.names", all_vertices)

# Activates/Deactivates vertices and edges based on your dataset's logic
for(i in 1:nrow(e0)) {
  if(e0$wt[i] == 1) {
    tail <- which(get.vertex.attribute(net_dynamic, "vertex.names") == e0$source[i])
    head <- which(get.vertex.attribute(net_dynamic, "vertex.names") == e0$target[i])
    if(length(tail) > 0 && length(head) > 0) {
      add.edge(net_dynamic, tail, head, onset=e0$time[i], terminus=e0$time[i]+1)
    }
  }
}

#Acctivates all vertices
activate.vertices(net_dynamic, onset = 0, terminus = length(wave_num) + 1)

#Get rid of the first 0-1 interval
deactivate.vertices(net_dynamic, onset = 0, length = 1)

#Create separate people not in the house for each waves
for(wave in wave_num){
  cur_wave <- paste("w", wave, sep = "")
  cur_name <- paste("notIn_w", wave, sep = "")
  individual_ids <- House_16_Presence %>%
    filter(!!sym(cur_wave) == 0) %>%
    pull(ID)
  assign(cur_name, individual_ids)
}

#Deactivate vertices who are not in the house
for(wave in wave_num){
  in_name <- paste("notIn_w", wave, sep = "")
  cur_active <- get(in_name)
  
  #Loop through the list of people not in the house
  for(i in 1:length(cur_active)){
    cur_id <- get_vertex_id_by_name(net_dynamic, cur_active[[i]])
    deactivate.vertices(net_dynamic, onset = wave, length = 1, v = cur_id, deactivate.edges = FALSE)
  }
}

# Active valid edges at each wave
for(wave in 1:length(wave_num)){
  #Find the current wave
  valid_edges_table_name <- paste("validEdge_w", wave, sep = "")
  valid_edges <- get(valid_edges_table_name)
  
  
  #Active edges at current wave
  for(i in 1:nrow(valid_edges)){
    source_vertex <- which(get.vertex.attribute(net_dynamic, "vertex.names") == valid_edges[i, "source"])
    target_vertex <- which(get.vertex.attribute(net_dynamic, "vertex.names") == valid_edges[i, "target"])
    
    #Check if the vertex exist
    if(length(source_vertex) > 0 && length(target_vertex) > 0){
      
      #Check if the edge already existed, add new edge and activatee if not
      existing_edge_ids <- get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex)
      if(length(existing_edge_ids) == 0){
        add.edge(net_dynamic, tail = source_vertex, head = target_vertex, activate = FALSE)
        activate.edges(net_dynamic, onset = wave, terminus = wave,
                       e = get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex))
      }
      
      #Activates an existing edge
      activate.edges(net_dynamic, onset = wave, terminus = wave,
                     e = get.edgeIDs(net_dynamic, v = source_vertex, alter = target_vertex))
    }
  }
}

# Assigns dynamic attributes (e.g., vertex color)
vertex_colors <- rep("cyan", network.size(net_dynamic))
names(vertex_colors) <- get.vertex.attribute(net_dynamic, "vertex.names")
vertex_colors["1228"] <- "red"  # Adjust based on actual vertex ID/names
set.vertex.attribute(net_dynamic, "vertex.color", vertex_colors)

# Renders the Dynamic Network
render.d3movie(net_dynamic,
               displaylabels = TRUE,
               vertex.col = "vertex.color",
               edge.lwd = 4,
               vertex.cex = 12,
               label = "vertex.names",
               output.mode = "html",
               filename = "dynamic_network_visualization.html")
