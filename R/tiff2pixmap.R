## tiff2pixmap
## read in TIFF image and output a Pixmap object.
## requires EBImage::readImage
## replaced deprecated rtiff:readTiff


tiff2pixmap=function(fn,pixmap=TRUE) {
        
        tiff=EBImage::readImage(fn)
        
        # d=imageData(tiff) # would need another import
        d=tiff@.Data
        
        w = dim(d)[1]
        h = dim(d)[2]
        

        # Replicating 2 dimensional matrix to create a 3 dimensional array
        if (length(dim(d))<3) d=replicate(3, d, simplify="array")
        
        # the way pixmap organizes the matrix appears to be transposed as
        # EBImage organzes it. to make it more pixmap way, need to transpose x,y
        # at each dimension
        dt=aperm(d, c(2,1,3))
        
        if (w > 0 && h > 0) {
            r= matrix(dt[,,1])
            g= matrix(dt[,,2])
            b= matrix(dt[,,3])
            
            rmx = max(r)
            gmx = max(g)
            bmx = max(b)

            
            if(pixmap) {
                pmap = pixmap::pixmapRGB(
                    data=array(data = c(r, g, b), dim = c(h, w, 3)), 
                    nrow=h, ncol=w,
                    bbox=NULL, bbcent=FALSE, cellres=c(1,1))
            } else {
                pmap = list(r = r, g=g, b=b)
            }

            
            return(pmap)
        } else {
            cat("Could not open", fn, ".  File corrupted or missing.\n")
        }
    }
