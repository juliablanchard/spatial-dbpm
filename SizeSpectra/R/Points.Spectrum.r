Points.Spectrum<-function(timestep, mass.lim, x.lim, y.lim,...){

  #Define entire mass spectrum
  mrange<-timestep@mrange
  xrange<-timestep@xrange
  yrange<-timestep@yrange

  #--------------------------------------------------------#
  # Produce three functions depending on Spatial Dimension #
  #--------------------------------------------------------#

  #---------------#
  # Aspatial Plot #
  #---------------#
  
  if(timestep@spatial_dim==0){
    if(missing(mass.lim)) mass.lim=vector() ; mass.lim[1] = min(mrange) ; mass.lim[2] = max(mrange)
    if(!is.numeric(mass.lim)) stop("Mass limits must be numeric")
    if(length(mass.lim)!=2) stop("Incorrect number of mass limits supplied")
    
    if(mass.lim[1]<min(mrange)) stop("Lower mass limit is too small")
    if(mass.lim[2]>max(mrange)) stop("Upper mass limit is too large")
    if(mass.lim[1]>=mass.lim[2]) stop("Lower mass limit must be smaller than Upper mass limit")
    if(mass.lim[1]==mass.lim[2]) stop("Mass limits must be distinct")
    if(!sum(mrange==mass.lim[1])) stop("Lower mass limit is not in discretisation")
    if(!sum(mrange==mass.lim[2])) stop("Upper mass limit is not in discretisation")
  
    #Finding index for minimum mass value to plot
    imin<-which(mrange==mass.lim[1])
    imax<-which(mrange==mass.lim[2])

    #Plot log plot
    uvals<-log(timestep@data[(imin:imax)+2])
    mvals<-mrange[imin:imax]
    points(mvals, uvals, ...)
  }

  #----------------#
  # 1D space Plots #
  #----------------#
  
  if(timestep@spatial_dim==1){
    if(missing(mass.lim)) {mass.lim=vector() ; mass.lim[1] = min(mrange) ; mass.lim[2] = max(mrange)}
    if(missing(x.lim)) {x.lim=vector() ; x.lim[1] = min(xrange) ; x.lim[2] = max(xrange)}
    if(!is.numeric(mass.lim)) stop("Mass limits must be numeric")
    if(!is.numeric(x.lim)) stop("X-Space limits must be numeric")
    if(length(mass.lim)>2) stop("Incorrect number of mass limits supplied")
    if(length(x.lim)>2) stop("Incorrect number of x-space limits supplied")
    
    if(length(mass.lim)==1 & length(x.lim)==1) stop("Nothing to plot") 
    if(length(mass.lim)==2 & length(x.lim)==2) stop("Cannot add to 3D plots")

    
    #Plot SizeSpectrum
    if(length(mass.lim)==2 & length(x.lim)==1){
      if(mass.lim[1]<min(mrange)) stop("Lower mass limit is too small")
      if(mass.lim[2]>max(mrange)) stop("Upper mass limit is too large")
      if(mass.lim[1]>=mass.lim[2]) stop("Lower mass limit must be smaller than upper mass limit")
      if(!sum(mrange==mass.lim[1])) stop("Lower mass limit is not in discretisation")
      if(!sum(mrange==mass.lim[2])) stop("Upper mass limit is not in discretisation")
      
      if(x.lim<min(xrange)) stop ("X-space value is too small")
      if(x.lim>max(xrange)) stop ("X-space value is too large")
      if(!sum(xrange==x.lim)) stop("X-space value is not in discretisation")

      #Finding index for minimum mass value to plot
      imin<-which(mrange==mass.lim[1])
      imax<-which(mrange==mass.lim[2])
      kval<-which(xrange==x.lim)
      
      uvals<-log(as.matrix(timestep@data)[kval,((imin:imax)+2)])
      mvals<-mrange[imin:imax]
      xvals<-xrange[kval]
      umin<-min(uvals[is.finite(uvals)])
      umax<-max(uvals[is.finite(uvals)])
      points(mvals, uvals,...)
    }
    
    #Plot Spatial Transect
    if(length(mass.lim)==1 & length(x.lim==2)){
      if(x.lim[1]<min(xrange)) stop("Lower x-space limit is too small")
      if(x.lim[2]>max(xrange)) stop("Upper x-space limit is too large")
      if(x.lim[1]>=x.lim[2]) stop("Lower x-space limit must be smaller than upper x-space limit")
      if(!sum(xrange==x.lim[1])) stop("Lower x-space limit is not in discretisation")
      if(!sum(xrange==x.lim[2])) stop("Upper x-space limit is not in discretisation")
      
      if(mass.lim<min(mrange)) stop ("Mass value is too small")
      if(mass.lim>max(mrange)) stop ("Mass value is too large")
      if(!sum(mrange==mass.lim)) stop("Mass value is not in discretisation")
      
      #Finding index for minimum mass value to plot
      ival<-which(mrange==mass.lim)
      kmin<-which(xrange==x.lim[1])
      kmax<-which(xrange==x.lim[2])
      
      uvals<-log(as.matrix(timestep@data)[(kmin:kmax),ival+2])
      mvals<-mrange[ival]
      xvals<-xrange[kmin:kmax]
      umin<-min(uvals[is.finite(uvals)])
      umax<-max(uvals[is.finite(uvals)])
      points(xvals, uvals,...)
    }
  }

  #----------------#
  # 2D space Plots #
  #----------------#  
  
  if(timestep@spatial_dim==2){
    if(missing(mass.lim)) {mass.lim=vector() ; mass.lim[1] = min(mrange) ; mass.lim[2] = max(mrange)}
    if(missing(x.lim)) {x.lim=vector() ; x.lim[1] = min(xrange) ; x.lim[2] = max(xrange)}
    if(missing(y.lim)) {y.lim=vector() ; y.lim[1] = min(yrange) ; y.lim[2] = max(yrange)}
    if(!is.numeric(mass.lim)) stop("Mass limits must be numeric")
    if(!is.numeric(x.lim)) stop("X-Space limits must be numeric")
    if(!is.numeric(y.lim)) stop("Y-Space limits must be numeric")
    if(length(mass.lim)>2) stop("Incorrect number of mass limits supplied")
    if(length(x.lim)>2) stop("Incorrect number of x-space limits supplied")
    if(length(y.lim)>2) stop("Incorrect number of y-space limits supplied")
    
    #Single point plot
    if(length(mass.lim)==1 & length(x.lim)==1 & length(y.lim)==1) stop("Nothing to plot")
    
    #Impossible plots
    if(length(mass.lim)==2 & length(x.lim)==2 & length(y.lim)==2) stop("Cannot plot 4D object")
    if(length(mass.lim)==2 & length(x.lim)==2 & length(y.lim)==1) stop("Cannot add to 3D plots")
    if(length(mass.lim)==2 & length(x.lim)==1 & length(y.lim)==2) stop("Cannot add to 3D plots")
    if(length(mass.lim)==1 & length(x.lim)==2 & length(y.lim)==2) stop("Cannot add to 3D plots")
    
    #Plot SizeSpectrum
    if(length(mass.lim)==2 & length(x.lim)==1 & length(y.lim)==1){
      if(mass.lim[1]<min(mrange)) stop("Lower mass limit is too small")
      if(mass.lim[2]>max(mrange)) stop("Upper mass limit is too large")
      if(mass.lim[1]>=mass.lim[2]) stop("Lower mass limit must be smaller than upper mass limit")
      if(!sum(mrange==mass.lim[1])) stop("Lower mass limit is not in discretisation")
      if(!sum(mrange==mass.lim[2])) stop("Upper mass limit is not in discretisation")
      
      if(x.lim<min(xrange)) stop ("X-space value is too small")
      if(x.lim>max(xrange)) stop ("X-space value is too large")
      if(!sum(xrange==x.lim)) stop("X-space value is not in discretisation")
      
      if(y.lim<min(yrange)) stop ("Y-space value is too small")
      if(y.lim>max(yrange)) stop ("Y-space value is too large")
      if(!sum(yrange==y.lim)) stop("Y-space value is not in discretisation")

      #Finding index for minimum mass value to plot
      imin<-which(mrange==mass.lim[1])
      imax<-which(mrange==mass.lim[2])
      kval<-which(xrange==x.lim)
      lval<-which(yrange==y.lim)
      
      uvals<-timestep@data[timestep@data$y==y.lim,]
      uvals<-log(as.matrix(uvals)[kval,((imin:imax)+2)])
      mvals<-mrange[imin:imax]
      xvals<-xrange[kval]
      yvals<-yrange[lval]
      points(mvals, uvals,...)
    }
    
    #Plot X-Spatial Transect
    if(length(mass.lim)==1 & length(x.lim==2) & length(y.lim)==1){
      if(x.lim[1]<min(xrange)) stop("Lower x-space limit is too small")
      if(x.lim[2]>max(xrange)) stop("Upper x-space limit is too large")
      if(x.lim[1]>=x.lim[2]) stop("Lower x-space limit must be smaller than upper x-space limit")
      if(!sum(xrange==x.lim[1])) stop("Lower x-space limit is not in discretisation")
      if(!sum(xrange==x.lim[2])) stop("Upper x-space limit is not in discretisation")
      
      if(mass.lim<min(mrange)) stop ("Mass value is too small")
      if(mass.lim>max(mrange)) stop ("Mass value is too large")
      if(!sum(mrange==mass.lim)) stop("Mass value is not in discretisation")

      if(y.lim<min(yrange)) stop ("Y-space value is too small")
      if(y.lim>max(yrange)) stop ("Y-space value is too large")
      if(!sum(yrange==y.lim)) stop("Y-space value is not in discretisation")
      
      #Finding index for minimum mass value to plot
      ival<-which(mrange==mass.lim)
      kmin<-which(xrange==x.lim[1])
      kmax<-which(xrange==x.lim[2])
      lval<-which(yrange==y.lim)
      
      uvals<-timestep@data[timestep@data$y==y.lim,]
      uvals<-log(as.matrix(uvals)[(kmin:kmax),ival+2])
      mvals<-mrange[ival]
      xvals<-xrange[kmin:kmax]
      yvals<-yrange[lval]
      points(xvals, uvals,...)
    }

    #Plot Y-Spatial Transect
    if(length(mass.lim)==1 & length(x.lim==1) & length(y.lim)==2){
      if(y.lim[1]<min(yrange)) stop("Lower y-space limit is too small")
      if(y.lim[2]>max(yrange)) stop("Upper y-space limit is too large")
      if(y.lim[1]>=y.lim[2]) stop("Lower y-space limit must be smaller than upper y-space limit")
      if(!sum(yrange==y.lim[1])) stop("Lower y-space limit is not in discretisation")
      if(!sum(yrange==y.lim[2])) stop("Upper y-space limit is not in discretisation")
      
      if(mass.lim<min(mrange)) stop ("Mass value is too small")
      if(mass.lim>max(mrange)) stop ("Mass value is too large")
      if(!sum(mrange==mass.lim)) stop("Mass value is not in discretisation")

      if(x.lim<min(xrange)) stop ("X-space value is too small")
      if(x.lim>max(xrange)) stop ("X-space value is too large")
      if(!sum(xrange==x.lim)) stop("X-space value is not in discretisation")
      
      #Finding index for minimum mass value to plot
      ival<-which(mrange==mass.lim)
      kval<-which(xrange==x.lim)
      lmin<-which(yrange==y.lim[1])
      lmax<-which(yrange==y.lim[2])
      
      uvals<-timestep@data[timestep@data$x==x.lim,]
      uvals<-log(as.matrix(uvals)[(lmin:lmax),ival+2])
      mvals<-mrange[ival]
      xvals<-xrange[kval]
      yvals<-yrange[lmin:lmax]
      points(xvals, uvals,...)
    }    
  }
}

