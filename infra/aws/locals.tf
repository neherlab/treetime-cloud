locals {
  environment  = "dev"
  cluster_name = "treetime-cloud-eks-dev-${random_string.suffix.result}"
}
