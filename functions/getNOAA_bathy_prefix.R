# custom version of getNOAA.bathy
# allows for file prefixs

getNOAA_bathy_prefix <- function (lon1, lon2, lat1, lat2, resolution = 4, keep = FALSE, 
					antimeridian = FALSE, prefix = "") 
{
	x1 = x2 = y1 = y2 = NULL
	if (lon1 < lon2) {
		x1 <- lon1
		x2 <- lon2
	}
	else {
		x2 <- lon1
		x1 <- lon2
	}
	if (lat1 < lat2) {
		y1 <- lat1
		y2 <- lat2
	}
	else {
		y2 <- lat1
		y1 <- lat2
	}
	res = resolution * 0.0166666666666667
	fetch <- function(x1, y1, x2, y2, res) {
		WEB.REQUEST <- paste("http://mapserver.ngdc.noaa.gov/cgi-bin/public/wcs/etopo1.xyz?filename=etopo1.xyz&request=getcoverage&version=1.0.0&service=wcs&coverage=etopo1&CRS=EPSG:4326&format=xyz&resx=", 
												 res, "&resy=", res, "&bbox=", x1, ",", y1, ",", x2, 
												 ",", y2, sep = "")
		print(WEB.REQUEST)
		dat <- suppressWarnings(try(read.table(WEB.REQUEST), 
																silent = TRUE))
		return(dat)
	}
	if (antimeridian) {
		FILE <- paste("marmap_coord_", x1, ";", y1, ";", x2, 
									";", y2, "_res_", resolution, "_anti", ".csv", sep = "")
	}
	else {
		FILE <- paste("marmap_coord_", x1, ";", y1, ";", x2, 
									";", y2, "_res_", resolution, ".csv", sep = "")
	}
	if (FILE %in% list.files(gsub("/", "", prefix))) {
		cat("File already exists ; loading '", FILE, "'", sep = "")
		exisiting.bathy <- read.bathy(paste0(prefix, FILE), header = T)
		return(exisiting.bathy)
	}
	else {
		if (antimeridian) {
			l1 <- x2
			l2 <- 180
			l3 <- -180
			l4 <- x1
			cat("Querying NOAA database ...\n")
			cat("This may take seconds to minutes, depending on grid size\n")
			left <- fetch(l1, y1, l2, y2, res)
			right <- fetch(l3, y1, l4, y2, res)
			if (is(left, "try-error") | is(right, "try-error")) {
				stop("The NOAA server cannot be reached\n")
			}
			else {
				cat("Building bathy matrix ...\n")
				left <- as.bathy(left)
				left <- left[-nrow(left), ]
				right <- as.bathy(right)
				rownames(right) <- as.numeric(rownames(right)) + 
					360
				bath2 <- rbind(left, right)
				class(bath2) <- "bathy"
				bath <- as.xyz(bath2)
			}
		}
		else {
			cat("Querying NOAA database ...\n")
			cat("This may take seconds to minutes, depending on grid size\n")
			bath <- fetch(x1, y1, x2, y2, res)
			if (is(bath, "try-error")) {
				stop("The NOAA server cannot be reached\n")
			}
			else {
				cat("Building bathy matrix ...\n")
				bath2 <- as.bathy(bath)
			}
		}
		if (keep) {
			write.table(bath, file = paste0(prefix, FILE), sep = ",", quote = FALSE, 
									row.names = FALSE)
		}
		return(bath2)
	}
}