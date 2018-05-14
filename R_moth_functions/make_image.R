make_image = function(FEM,habitat,kval)  
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
  # vectors X and Y;
  #
  
  #   if (!is.fd(fdobj))
  #   {
  #     stop('FDOBJ is not an FD object')
  #   }
  
  nodes = FEM$FEMbasis$mesh$nodes
  triangles = FEM$FEMbasis$mesh$triangles
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    #rgl.open()
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    nhab<-length(kval)
    for (hab in 1:nhab){
      sel=habitat==hab
      z = kval[hab]*coeff
      if (hab==1){
        miz=min(z)
        maz=max(z)
      } else
      {miz=min(miz,min(z))
      maz=max(maz,max(z))
      }
    }
    for (hab in 1:nhab){
      sel=habitat==hab
      z = kval[hab]*coeff[as.vector(t(triangles[sel,])),isurf]
      rgl.triangles(x = nodes[as.vector(t(triangles[sel,])) ,1], y = nodes[as.vector(t(triangles[sel,])) ,2], 
                    z=0, 
                    color = heat[round(99*(z- miz)/(maz-miz))+1])
    }
    aspect3d(2,2,1)
    rgl.viewpoint(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}
