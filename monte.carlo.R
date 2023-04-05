simulate = function(adj.mat, p = 0.2, n.iter = 5000){
    if (p <= 0 | p >= 1){
        stop("Invalid limits. 0 < limit < 1 is valid")
    }
    if (p >= 0.5){
        p.mark = 1 - p
    }
    if (p < 0.5){
        p.mark = p
    }
    histogram = data.frame(Size = 0,SumProbabilities = 0, MeanProbabilities = 0,ExpectedNn = 0,PMaxNGreaterNi = 0)
    max.clust.size = data.frame()
    pb <- progress::progress_bar$new(
        format = "process [:bar] :current in :elapsed :eta",
        total = n.iter, clear = FALSE, width= 60)
    for (i in 1:n.iter) {
        #mark random regions with probability marks.probs
        itendif = sample(0:1,nrow(adj.mat),replace = T,prob = c(1 - p.mark, p.mark))
        select = which(itendif != 0, arr.ind = T)
        #find connected(neibourhood AND marked) regions - clusters and this sizes
        graph = igraph::graph.adjacency(adj.mat[select, select])
        comps <- igraph::components(graph)
        clusters = igraph::groups(comps)
        sizes = sapply(clusters, length)

        if (length(sizes) == 0) {
            max.clust.size[nrow(max.clust.size) + 1, 1] = 0
            size.probs = 0
            names(size.probs) = c("0")
        } 
        else{
            max.clust.size[nrow(max.clust.size)+1,1] = max(sizes)
            num.sizes = table(sizes)
            size.probs = num.sizes / length(sizes)
        }
            
        size.probs = size.probs[order(as.numeric(names(size.probs)))]
        k = 1
        for (j in names(size.probs)) {
            id = as.numeric(j)
            if(!(id %in% histogram$Size)){
                histogram[nrow(histogram)+1,] = c(id,0,0,num.sizes[j],0)
            }
                
            histogram[histogram$Size == id, "SumProbabilities"] = histogram[histogram$Size == id, "SumProbabilities"] + size.probs[j]
            histogram[histogram$Size == id, "ExpectedNn"] = histogram[histogram$Size == id, "ExpectedNn"] + num.sizes[j]
            revCDF = ifelse(length(size.probs) == 1, size.probs, sum(size.probs[k:nrow(size.probs)]))
            k = k + 1
            histogram[histogram$Size == id, "PMaxNGreaterNi"] = histogram[histogram$Size == id, "PMaxNGreaterNi"] + revCDF
        }
        pb$tick()
    }
    #PDF of cluster size
    histogram$MeanProbabilities = histogram$SumProbabilities / n.iter
    histogram$ExpectedNn = histogram$ExpectedNn / n.iter
    histogram$PMaxNGreaterNi = histogram$PMaxNGreaterNi / n.iter
    histogram = histogram[order(histogram$Size), ]
    return(histogram[-1,])
}

plot.sim = function(hist){
    `%>%` <- magrittr::`%>%`
    fig1 = plotly::plot_ly(hist, x = ~Size, y = ~MeanProbabilities,type = 'bar',color = I("blue"), alpha = 0.85) %>%
        plotly::layout(xaxis = list(title="Size"),yaxis = list(title="Probability"), showlegend = F)
    fig2 = plotly::plot_ly(hist, x = ~Size, y = ~PMaxNGreaterNi,type = 'bar',color = I("blue"), alpha = 0.85) %>%
        plotly::layout(xaxis = list(title="Size"),yaxis = list(title="Probability"), showlegend = F)
    annotations = list( 
        list( 
            x = 0.2,  
            y = 1.0,  
            text = "Size distribution",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "bottom",  
            yshift = -10,
            showarrow = FALSE 
        ),  
        list( 
            x = 0.8,  
            y = 1,  
            text = "Maximal size distribution",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "bottom", 
            yshift = -10,
            showarrow = FALSE 
        ))
    fig = plotly::subplot(fig1,fig2,shareX = T, shareY = T) %>% plotly::layout(annotations = annotations)
    fig
}

build.critical.region = function(simulation, alpha = 0.05){
    crit = simulation
    a1 = alpha / 2
    a2 = alpha / 2
    crit$distToCrit = abs(crit$PMaxNGreaterNi - a2)
    ind = which.min(crit$distToCrit)
    a2 = crit$PMaxNGreaterNi[ind]
    a1 = alpha - a2
    n.crit = crit$Size[ind]
    a.uni = a1 / (n.crit - 1)
    f.cr = function(s.cr,size,e.n){return(pchisq(s.cr,df = 2 * size) - (1 - a.uni / e.n))}
    tol = 0.00001
    
    for (j in 1:nrow(crit)) {
        if(crit[j,"Size"] <= n.crit){
            crit[j,"S.crit"] = rootSolve::uniroot.all(f.cr,
                                                        c(qchisq(tol,2*crit[j,"Size"]), qchisq(1-tol,2*crit[j,"Size"])),
                                                        tol = tol,
                                                        size = crit[j,"Size"],
                                                        e.n = crit[j,"ExpectedNn"])[1]
        }
        else{
            crit[j,"S.crit"] = 0
        }
    }
    crit = rbind.data.frame(crit,data.frame(Size = n.crit,
                                                SumProbabilities = 0,
                                                MeanProbabilities = 0,
                                                ExpectedNn = 0,
                                                PMaxNGreaterNi = 0,
                                                distToCrit = 0,
                                                S.crit = 0))
    crit = crit[order(crit$Size),]
    return(crit)
}

plot.crit.region = function(crit){
    `%>%` <- magrittr::`%>%`
    fig = plotly::plot_ly(crit, x = ~Size, y = ~S.crit,type = 'scatter', mode = 'lines+markers',color = I("black"), name = "Критическая область") %>%
        plotly::layout(xaxis = list(title="Size"),yaxis = list(title="Statistic"))
    fig
}

search.structures = function(map, adj.mat, data, column, crit, simulation, key = "GID_1", clusters = T, p.limit = 0.8){
    map.data = dplyr::inner_join(map@data, data, by = key)
    if (clusters){
        s.n = function(x, p, data){
            return(-2 * sum((log((1 - data[which(data[, key] %in% x),column]+0.0001) / (1 - p)))))
        }
        map.data$CRIT = ifelse(map.data[, column] >= p.limit, 1, 0)
    }
    else{
        s.n.down = function(x, p, data){
            return(-2 * sum(log((data[which(data[, key] %in% x),column]+0.0001) / p)))
        }
        map.data$CRIT = ifelse(map.data[, column] <= p.limit, 1, 0)
    }
    select = which(map.data$CRIT == 1,arr.ind = T)
    if(length(select) != 0){
        graph <- igraph::graph.adjacency(adj.mat[select,select])
        comps <- igraph::components(graph)
        clusters = igraph::groups(comps)
        if(length(clusters) > 0){
            s = sapply(clusters,s.n,p = p.limit, data = map.data)
            sizes = sapply(clusters, length)
            cls = as.data.frame(table(sizes))
            structures = data.frame(ID = names(s),N = sizes,S = s)
            #names(clusters) = data.discharges$ID
            for (j in 1:nrow(structures)) {
                structures$Label[j] = paste(clusters[[as.character(structures$ID[j])]],collapse = "; ")
                structures$Names[j] = paste(map.data$NAME_1[which(map.data$GID_1 %in% clusters[[as.character(structures$ID[j])]])],collapse = "; ")
            }     
            nearest = unlist(lapply(structures$N,function(x){which.min(abs(crit$Size - x))}))
            structures$S.crit = crit[nearest,"S.crit"]
            structures$N.crit = crit$Size[which.min(crit$distToCrit)] 
            return(structures)
        }
        else{
            return(NULL)
        }
    }
    else{
        return(NULL)
    }
}

get.significant.structures = function(structures){
    `%>%` <- magrittr::`%>%`
    structures %>% 
        dplyr::filter(N >= N.crit | S >= S.crit)
}

plot.result = function(crit, structures, significant = NULL){
    crit.region = plot.crit.region(crit)
    `%>%` <- magrittr::`%>%`
    structures.points = crit.region %>%
        plotly::add_trace(data = structures,x = ~N,y = ~S, name = 'Незначимые', mode = 'markers',color = I("blue"), size = 3, text = ~Label)
    if(!is.null(significant)){
        structures.points = structures.points %>% 
            plotly::add_trace(data = significant,x = ~N,y = ~S, name = 'Значимые', mode = 'markers',color = I("black"), size = 3, text = ~Label)
    }
    structures.points
}  

plot.map.structures = function(map, structures = NULL, data, key = "GID_1",render.column,range = c(0,1)){
    map1 = sf::st_shift_longitude(sf::st_as_sf(map))
    map1 = as(map1, "Spatial")
    map1@data = dplyr::inner_join(map1@data, data, by = key)
    seq = seq(from = range[1], to = range[2], 
              by = (range[2] - range[1]) / 1000)
    p <- leaflet::colorNumeric(grDevices::colorRampPalette(c("green", "red"))(length(seq)),domain = seq)
    coords = sp::coordinates(map1)
    l = leaflet::leaflet() %>% 
        leaflet::addProviderTiles("CartoDB.Positron") %>% 
        leaflet::addPolygons(
            data = map1,
            weight = 1,
            fillOpacity = 0.25,
            color = p(map1@data[,render.column]),
            layerId = map1@data$GID_1,
            highlight = leaflet::highlightOptions(
                weight = 5,
                color = p(map1@data[,render.column]),
                fillOpacity = 0.7,
                bringToFront = TRUE
            )
        ) %>% 
        leaflet::addLegend(pal = p,
                  title = render.column,
                  values = map1@data[,render.column],
                  layerId = "Legend") #%>% 
        #leaflet::fitBounds(min(coords[,1]),min(coords[,2]),max(coords[,1]),max(coords[,2]))
    if(!is.null(structures)){
        str.ids = stringr::str_split(structures$Label, "; ")
        names(str.ids) =  structures$ID
        for (i in structures$ID) {
            if(structures[structures$ID == i,"N"] >= structures[structures$ID == i,"N.crit"] | 
               structures[structures$ID == i,"S"] >= structures[structures$ID == i,"S.crit"]){
                color = "black"
            }
            else{
                color = "blue"
            }
            region = map1[which(map1@data$GID_1 %in% str.ids[[i]]),]
            #print(which(copy.map()@data$GID_1 %in% clusters[[i]]))
            union = surveillance::unionSpatialPolygons(region)
            l = l %>%
                leaflet::addPolylines(
                    data = union,
                    group = "union",
                    color = color,
                    weight = 5,
                    opacity = 1
                )
        }
    }
    l
}

calc.spatial.standard = function(data){
    M = mean(data,na.rm = T)
    X = (data - M) / sqrt(M)
    pnorm(X)
}

calc.spatial.fisher = function(data, population){
    f = data / population
    f.average = sum(data, na.rm = T) / sum(population, na.rm = T)
    z = 2 * sqrt(population) * asin(f - f.average)
    pnorm(z)
}

calc.time.standard = function(t1, t2){
    X = (t2 - t1) / (sqrt(t1 + t2))
    pnorm(X)
}
        
library(dplyr)
map = sf::st_read("~/Documents/Pyat/test_art5/HIV_RUS", layer = "Russia")
map = as(map,"Spatial")
map = map[(map@data$GID_1 != "RUS.84_1") &(map@data$GID_1 != "RUS.85_1"),]
adjacencyMatrix = surveillance::poly2adjmat(map, queen = T, zero.policy = T,row.names = map@data$GID_1)
new.map = rgdal::readOGR("~/Documents/Pyat/test_art5/New_RUS/", layer = 'Russia')
map@data = inner_join(map@data, new.map@data %>% dplyr::select(contains("pop"),"GID_1"),by = "GID_1")
data = map@data %>% dplyr::select(all_of(c("GID_1", "NAME_1", "Y2011", "Y2012", "popY2011", "popY2012"))) %>% 
    mutate(Y2011 = Y2011 / 100000 * popY2011,
           Y2012 = Y2012 / 100000 * popY2012)
data = data %>% 
    mutate(p.2011 = calc.spatial.fisher(Y2011, popY2011),
           p.2012 = calc.spatial.fisher(Y2012, popY2012),
           p.delta = calc.time.standard(Y2011 / popY2011 * 100000,
                                        Y2012 / popY2012 * 100000))
    #        p.delta = pnorm(delta)) %>% 
    # dplyr::select(all_of(c("GID_1", "p.2011","p.2012","p.delta")))

up = simulate(adjacencyMatrix, p = 0.8)
crit.up = build.critical.region(up)
str.up = search.structures(map, adjacencyMatrix, data, "p.delta",crit.up, up,  p.limit = 0.8)
signif.up = get.significant.structures(str.up)
plot.result(crit.up,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.delta")

map = sf::st_read("~/Documents/Pyat/test_art5/HIV_RUS", layer = "Russia")
map = as(map,"Spatial")
map = map[(map@data$GID_1 != "RUS.84_1") &(map@data$GID_1 != "RUS.85_1"),]
adjacencyMatrix = surveillance::poly2adjmat(map, queen = T, zero.policy = T,row.names = map@data$GID_1)
new.map = rgdal::readOGR("~/Documents/Pyat/test_art5/New_RUS/", layer = 'Russia')
map@data = inner_join(map@data, new.map@data %>% dplyr::select(contains("pop"),"GID_1"),by = "GID_1")

data = map@data %>% 
    mutate(Y2011_new = Y2011 / 100000 * popY2011,
           Y2012_new = Y2012 / 100000 * popY2012,
           Y2013_new = Y2013 / 100000 * popY2013,
           p.2011 = calc.spatial.fisher(Y2011_new, popY2011),
           p.2012 = calc.spatial.fisher(Y2012_new, popY2012),
           p.2013 = calc.spatial.fisher(Y2013_new, popY2013),
           p.delta = calc.time.standard(Y2011, Y2012),
           p.delta1 = calc.time.standard(Y2012, Y2013))

up = simulate(adjacencyMatrix, p = 0.85)
crit.up = build.critical.region(up)

str.up = search.structures(map, adjacencyMatrix, data, "p.2011",crit.up, up,  p.limit = 0.85)
signif.up = get.significant.structures(str.up)
plot.result(crit.up,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.2011")

str.up = search.structures(map, adjacencyMatrix, data, "p.2012",crit.up, up,  p.limit = 0.85)
signif.up = get.significant.structures(str.up)
plot.result(crit.up,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.2012")

str.up = search.structures(map, adjacencyMatrix, data, "p.2013",crit.up, up,  p.limit = 0.85)
signif.up = get.significant.structures(str.up)
plot.result(crit.up,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.2013")


up1 = simulate(adjacencyMatrix, p = 0.85)
crit.up1 = build.critical.region(up1)
str.up = search.structures(map, adjacencyMatrix, data, "p.delta",crit.up1, up1,  p.limit = 0.85)
signif.up = get.significant.structures(str.up)
plot.result(crit.up1,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.delta")


str.up = search.structures(map, adjacencyMatrix, data, "p.delta1",crit.up1, up1,  p.limit = 0.85)
signif.up = get.significant.structures(str.up)
plot.result(crit.up1,str.up,signif.up)
plot.map.structures(map, str.up,data = data, render.column = "p.delta1")
