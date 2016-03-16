
#
# functions for plotting a tree
#

# modified by kieran samuk mar 2016
# note: these functions are rather inscrutable

library(RColorBrewer)
set_y_coords = function(d){
	i = which(d[,3]=="ROOT")
	y = d[i,8]/ (d[i,8]+d[i,10])
	d[i,]$y = 1-y
	d[i,]$ymin = 0
	d[i,]$ymax = 1
	c1 = d[i,7]
	c2 = d[i,9]
	ni = which(d[,1]==c1)
	ny = d[ni,8]/ (d[ni,8]+d[ni,10])
        d[ni,]$ymin = 1-y
	d[ni,]$ymax = 1
	d[ni,]$y = 1- ny*(y)

        ni = which(d[,1]==c2)
        ny = d[ni,8]/ (d[ni,8]+d[ni,10])
        d[ni,]$ymin = 0
        d[ni,]$ymax = 1-y
        d[ni,]$y = (1-y)-ny*(1-y)

	for (j in 1:nrow(d)){
		d = set_y_coord(d, j)
	}	
	return(d)
}

set_y_coord = function(d, i){
	index = d[i,1]
	parent = d[i,6]
	if (!is.na(d[i,]$y)){
		return(d)
	}
	tmp = d[d[,1] == parent,]
	if ( is.na(tmp[1,]$y)){
		d = set_y_coord(d, which(d[,1]==parent))
		tmp = d[d[,1]== parent,]
	}
	py = tmp[1,]$y
	pymin = tmp[1,]$ymin
	pymax = tmp[1,]$ymax
	f = d[i,8]/( d[i,8]+d[i,10])
	#print (paste(i, index, py, pymin, pymax, f))
	if (tmp[1,7] == index){
		d[i,]$ymin = py
		d[i,]$ymax = pymax
		d[i,]$y = pymax-f*(pymax-py)
		if (d[i,5]== "TIP"){
			d[i,]$y = (py+pymax)/2
		}
	}
	else{
		d[i,]$ymin = pymin
		d[i,]$ymax = py
		d[i,]$y = py-f*(py-pymin)
                if (d[i,5]== "TIP"){
                        d[i,]$y = (pymin+py)/2
                }	
	
	}
	return(d)
}


set_x_coords = function(d, e){
	i = which(d[,3]=="ROOT")
	index = d[i,1]
        d[i,]$x = 0
        c1 = d[i,7]
        c2 = d[i,9]
        ni = which(d[,1]==c1)
	tmpx =  e[e[,1]==index & e[,2] == c1,3]
	if (length(tmpx) == 0){
		tmp = e[e[,1] == index,]
		tmpc1 = tmp[1,2]
		if ( d[d[,1]==tmpc1,4] != "MIG"){
			tmpc1 = tmp[2,2]
		}
		tmpx = get_dist_to_nmig(d, e, index, tmpc1)
	}
	if(tmpx < 0){
		tmpx = 0
	}
        d[ni,]$x = tmpx

        ni = which(d[,1]==c2)
        tmpx =  e[e[,1]==index & e[,2] == c2,3]
	if (length(tmpx) == 0){
    		tmp = e[e[,1] == index,]
                tmpc2 = tmp[2,2]
                if ( d[d[,1]==tmpc2,4] != "MIG"){
                        tmpc2 = tmp[1,2]
                }
                tmpx = get_dist_to_nmig(d, e, index, tmpc2)
	}
        if(tmpx < 0){
                tmpx = 0
        }
        d[ni,]$x = tmpx
        
        for (j in 1:nrow(d)){
                d = set_x_coord(d, e, j)
        }
        return(d)
	#print(d)
}


set_x_coord = function(d, e, i){
	index = d[i,1]
        parent = d[i,6]
        if (!is.na(d[i,]$x)){
                return(d)
        }
        tmp = d[d[,1] == parent,]
        if ( is.na(tmp[1,]$x)){
                d = set_x_coord(d, e, which(d[,1]==parent))
                tmp = d[d[,1]== parent,]
        }
	#print (paste(parent, index))
	tmpx =  e[e[,1]==parent & e[,2] == index,3]
	if (length(tmpx) == 0){
                tmp2 = e[e[,1] == parent,]
                tmpc2 = tmp2[2,2]
		#print
                if ( d[d[,1]==tmpc2,4] != "MIG"){
                        tmpc2 = tmp2[1,2]
                }
		tmpx = get_dist_to_nmig(d, e, parent, tmpc2)
	}
        if(tmpx < 0){
                tmpx = 0
        }
	d[i,]$x = tmp[1,]$x+ tmpx
	return(d)
}



set_mig_coords = function(d, e){
	for (j in 1:nrow(d)){
		if (d[j,4] == "MIG"){
			p = d[d[,1] == d[j,6],]
			c = d[d[,1] == d[j,7],]
			tmpe = e[e[,1] == d[j,1],]
			y1 = p[1,]$y
			y2 = c[1,]$y
			x1 = p[1,]$x
			x2 = c[1,]$x

			mf = tmpe[1,6]	
			if (is.nan(mf)){
				mf = 0
			}
			#d[j,]$y = (y1+y2)* mf
                        #d[j,]$x = (x1+x2) *mf
                        d[j,]$y = y1+(y2-y1)* mf
			#print(paste(mf, x1, x2))
                        d[j,]$x = x1+(x2-x1) *mf
		}	

	}
	return(d)

}

get_f = function(stem){
	d = paste(stem, ".cov.gz", sep = "")
	d2 = paste(stem, ".modelcov.gz", sep = "")
	d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
	d2 = read.table(gzfile(d2), as.is = T, comment.char = "", quote = "")
	d = d[order(names(d)), order(names(d))]
	d2 = d2[order(names(d2)), order(names(d2))]
	tmpcf = vector()
        tmpmcf = vector()
        for (j in 1:nrow(d)){
                for (k in (j+1):nrow(d)){
                        tmpcf = append(tmpcf, d[j,k])
                        tmpmcf = append(tmpmcf, d[j,k] - d2[j,k])
                }
        }
        tmpv = var(tmpmcf)/var(tmpcf)
	return(1-tmpv)

}



get_dist_to_nmig = function(d, e, n1, n2){
	toreturn = e[e[,1] == n1 & e[,2] == n2,3]
	#print(toreturn)
	while ( d[d[,1] ==n2,4] == "MIG"){
		tmp = e[e[,1] == n2 & e[,5] == "NOT_MIG",]
		toreturn = toreturn+tmp[1,3]
		n2 = tmp[1,2]
	}
	return(toreturn)
}

flip_node = function(d, n){
	i = which(d[,1] == n)
	t1 = d[i,7]
	t2 = d[i,8]
	d[i,7] = d[i,9]
	d[i,8] = d[i,10]
	d[i,9] = t1
	d[i,10] = t2
	return(d)

}

# wrote this so I can inject arbitrary info into the d object before it gets passed to treemix internal
get_treemix_d_object = function(stem, o = NA, flip = vector()){
	d = paste(stem, ".vertices.gz", sep = "")
	e = paste(stem, ".edges.gz", sep = "")
	se = paste(stem, ".covse.gz", sep = "")
	d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
	e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
	if (!is.na(o)){
		o = read.table(o, as.is = T, comment.char = "", quote = "")
	}
	e[,3] = e[,3]*e[,4]
	e[,3] = e[,3]*e[,4]
	
	se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
	m1 = apply(se, 1, mean)
	m = mean(m1)
	#m = 0
	for(i in 1:length(flip)){
		d = flip_node(d, flip[i])
	}
	d$x = "NA"
	d$y = "NA"
	d$ymin = "NA"
	d$ymax = "NA"
	d$x = as.numeric(d$x)
	d$y = as.numeric(d$y)
	d$ymin = as.numeric(d$ymin)
	d$ymax = as.numeric(d$ymax)

	d = set_y_coords(d)
	d = set_x_coords(d, e)
	#print(d)
	d = set_mig_coords(d, e)
	
	return(d)
}

get_dist_to_nmig = function(d, e, n1, n2){
	toreturn = e[e[,1] == n1 & e[,2] == n2,3]
	#print(toreturn)
	while ( d[d[,1] ==n2,4] == "MIG"){
		tmp = e[e[,1] == n2 & e[,5] == "NOT_MIG",]
		toreturn = toreturn+tmp[1,3]
		n2 = tmp[1,2]
	}
	return(toreturn)
}

flip_node = function(d, n){
	i = which(d[,1] == n)
	t1 = d[i,7]
	t2 = d[i,8]
	d[i,7] = d[i,9]
	d[i,8] = d[i,10]
	d[i,9] = t1
	d[i,10] = t2
	return(d)

}

plot_modelcov = function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T){
        c = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T)
        o = read.table(pop_order, as.is = T, comment.char = "", quote = "")


        toplot = data.frame(matrix(nrow = nrow(c), ncol = ncol(c)))
        for(i in 1:nrow(o)){
                for( j in 1:nrow(o)){

                        toplot[i, j] = c[which(names(c)==o[i,1]), which(names(c)==o[j,1])]
                }
        }
        if (usemax){
                m1 = max(abs(toplot))
                max = m1*1.1
                min = -(m1*1.1)
        }
        names(toplot) = o[,1]
        plot_resid_internal(toplot, max = max, min = min)
}



plot_cov = function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T, wcols = ""){
        c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T)
        o = read.table(pop_order, as.is = T)


        toplot = data.frame(matrix(nrow = nrow(c), ncol = ncol(c)))
        for(i in 1:nrow(o)){
                for( j in 1:nrow(o)){

                        toplot[i, j] = c[which(names(c)==o[i,1]), which(names(c)==o[j,1])]
                }
        }
        if (usemax){
                m1 = max(abs(toplot))
                max = m1*1.1
                min = 0
        }
        names(toplot) = o[,1]
        plot_cov_internal(toplot, max = max, min = min, wcols = wcols, o = o, cex = cex)
}


plot_resid = function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T, wcols = "r"){
	c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
	m = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
	names(c) = rownames(c)
	names(m) = rownames(m)
	o = read.table(pop_order, as.is = T, comment.char = "", quote = "")
	se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
	mse = apply(se, 1, mean)
	mse = mean(mse)
	print(mse)	
	c = c[order(names(c)), order(names(c))]
	m = m[order(names(m)), order(names(m))]
	tmp = c -m 
	#tmp = m - c
	#tmp = (m-c)/m
	#print(tmp)
	toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
	for(i in 1:nrow(o)){
	        for( j in 1:nrow(o)){
			#print(paste(o[i,1], o[j,1]))
			if (o[i,1] %in% names(tmp) ==F){
				print(paste("not found", o[i,1]))
			}
			if (o[j,1] %in% names(tmp) ==F){
				print(paste("not found", o[j,1]))
			}
        	        toplot[i, j] = tmp[which(names(tmp)==o[i,1]), which(names(tmp)==o[j,1])]
        	}
	}
	#print(toplot)
	if (usemax){
		m1 = max(abs(toplot), na.rm = T)
		max = m1*1.02
		min = -(m1*1.02)	
	}
	print("here")
	names(toplot) = o[,1]
	toreturn = plot_resid_internal(toplot, max = max, min = min, wcols = wcols, mse = mse, o = o, cex = cex)
	return(toreturn)
}

plot_cov_internal = function(d, o = NA, max = 0.009, min = -0.009, cex =0.5, wcols = "", mse = 5){
        npop = nrow(d)
        width = 1/npop
        height = 1/npop
        colors = brewer.pal(9, "Spectral")
        colors = c("red", "orange","yellow", "white", "green", "blue", "black")
        pal = colorRampPalette(colors)
        ncol = 80
        cols = pal(ncol)
        plot("NA", xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = "", ylab = "")
        for (i in 1:npop){
                for( j in 1:i){
                        v = d[i,j]
                        col= "white"
                        if (v < 0){
                                if (wcols == "rb"){
                                col = rgb(0, 0, 1, v/min)
                                }
                                else{
                                #col = rgb(0, 0, 1, 0.1+0.9*(v/min))
                                col = cols[ncol/2-floor( (v/min)*(ncol/2))]
                                }
                        }
                        else{
                                if (wcols == "rb"){
                                col = rgb(1, 0, 0, v/max)
                                }
                                else{
                                #col = rgb(1, 0, 0, 0.1+0.9*(v/max))
                                col = cols[ceiling((v/max)*(ncol))]
                                }
                        }
                        xmin = j/npop - 1/npop
                        xmax = j/npop
                        ymin = 1-(i/npop)
                        ymax = 1-(i/npop)+1/npop
			if (v == 0){ col = "white"}
                        rect(xmin, ymin, xmax, ymax, col = col, border = col)
                }
                tcol = "black"
                tmp = o[o[,1] == names(d)[i],]
                if (length(tmp) != 1){
                        tcol = tmp[1,2]
                }
                mtext(names(d)[i], side = 2, at = 1-i/npop+0.5/npop, las = 1, cex = cex, col = tcol)
                mtext(names(d)[i], side = 1, at =  i/npop-0.5/npop, las = 3, cex = cex, col = tcol)
        }
        if ( !is.na(mse)){
                ymi = 0.5
                yma = 0.9
                w = (yma-ymi)/ncol
                xma = 0.80
                lmi = round(min, digits = 1)
                lma = round(max, digits = 1)
                print(cols)
                print(ymi+(0:ncol)*w)
                rect( rep(0.75, ncol), ymi+(0:(ncol-1))*w, rep(xma, ncol), ymi+(1:ncol)*w, col = cols, border = cols)
                text(xma+0.01, ymi, lab = paste(lmi),  adj = 0, cex = 0.8)
                text(xma+0.01, yma, lab = paste(lma, "(Variance)"), adj = 0, cex = 0.8)

        }
        return(d)
        #image(as.matrix(d), col = cols)
}

plot_resid_internal = function(d, o = NA, max = 0.009, min = -0.009, cex =0.5, wcols = "rb", mse = NA){
        npop = nrow(d)
        width = 1/npop
        height = 1/npop
	colors = brewer.pal(9, "Spectral")
	colors = c("red", "orange","yellow", "white", "green", "blue", "black")
	pal = colorRampPalette(colors)
	ncol = 80
	cols = pal(ncol)
        plot("NA", xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = "", ylab = "")
        for (i in 1:npop){
                for( j in 1:i){
                        v = d[i,j]
			print(paste(i, j, v))
                        col= "white"
                        if (v < 0){
				if (wcols == "rb"){
				col = rgb(0, 0, 1, v/min)
				}
				else{
                                #col = rgb(0, 0, 1, 0.1+0.9*(v/min))
				col = cols[ncol/2-floor( (v/min)*(ncol/2))]
				#col = "white"
				}
                        }
                        else{
				if (wcols == "rb"){
				col = rgb(1, 0, 0, v/max)
				}
				else{
                                #col = rgb(1, 0, 0, 0.1+0.9*(v/max))
				col = cols[ncol/2+ceiling((v/max)*(ncol/2))]
				}
                        }
                        xmin = j/npop - 1/npop
                        xmax = j/npop
                        ymin = 1-(i/npop)
                        ymax = 1-(i/npop)+1/npop
                        rect(xmin, ymin, xmax, ymax, col = col, border = col)
                }
		tcol = "black"
		tmp = o[o[,1] == names(d)[i],]
		if (length(tmp) != 1){
			tcol = tmp[1,2]
		}
                mtext(names(d)[i], side = 2, at = 1-i/npop+0.5/npop, las = 1, cex = cex, col = tcol)
                mtext(names(d)[i], side = 1, at =  i/npop-0.5/npop, las = 3, cex = cex, col = tcol)
        }
	if ( !is.na(mse)){
                ymi = 0.5
                yma = 0.9
                w = (yma-ymi)/ncol
                xma = 0.80
		lmi = round(min/mse, digits = 1)
		lma = round(max/mse, digits = 1)
		print(cols)
		print(ymi+(0:ncol)*w)
                rect( rep(0.75, ncol), ymi+(0:(ncol-1))*w, rep(xma, ncol), ymi+(1:ncol)*w, col = cols, border = cols)
                text(xma+0.01, ymi, lab = paste(lmi, "SE"),  adj = 0, cex = 0.8)
                text(xma+0.01, yma, lab = paste(lma, "SE"), adj = 0, cex = 0.8)

	}
	return(d)
	#image(as.matrix(d), col = cols)
}


