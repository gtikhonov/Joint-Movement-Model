make_FEM_basis = function(file)
{
  load(file=file)
  mesh = create.MESH.2D(nodes=my.tp$P,segments=my.tp$S,triangles=my.tp$T)
  habitat = my.tp$VA[,1]+1
  MRR = my.tp$VA[,2]
  #plot(mesh)
  FEMbasis = create.FEM.basis(mesh)
  nnodes = nrow(mesh$nodes)
  nMRR = max(MRR)
  MRRM = matrix(0,nnodes,nMRR) #Matrix showing the nodes that belong to each MRR site
  for (StartingTrap in 1:nMRR){
    MRRM[unique(mesh$triangles[MRR==StartingTrap]),StartingTrap] = 1  
  }
  multiple = which(rowSums(MRRM)>1)
  for (i in 1:length(multiple)){
    cases = which(MRRM[multiple[i],]==1)
    new = rep(0,nMRR)
    new[sample(cases,1)] = 1
    MRRM[multiple[i],] = new
  }
  triangles = FEMbasis$mesh$triangles
  nele  = nrow(triangles)
  nele_local = ncol(triangles)^2
  Agrid = matrix(0,nele*nele_local,2) 
  for (el in 1:nele){
    ind = triangles[el,]
    Agrid[nele_local*(el-1)+(1:nele_local),] = as.matrix(expand.grid(ind,ind))
  }
  FEMbasis$Agrid = Agrid
  FEMbasis$habitat = habitat
  FEMbasis$MRRM = Matrix(MRRM)
  return(FEMbasis)
}
