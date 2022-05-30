library(EthSEQ)

ethseq.Analysis(
  bam.list = "./BAMs_List.txt",
  out.dir = "./",
  model.gds = "./SS2.Light.Model.gds",
  cores=1,
  mbq = 20,
  mrq = 1,
  mdc = 10,
  run.genotype = TRUE,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space="3D")
