module "eks" {
  source  = "terraform-aws-modules/eks/aws"
  version = "13.2.1"

  cluster_name       = local.cluster_name
  cluster_version    = "1.18"
  subnets            = module.vpc.private_subnets
  vpc_id             = module.vpc.vpc_id
  manage_aws_auth    = true
  config_output_path = "kubeconfig/kubeconfig.yaml"

  worker_groups = [
    {
      name                          = "worker-group-basic"
      instance_type                 = "t2.small"
      spot_price                    = ""
      asg_min_size                  = 1
      asg_max_size                  = 1
      asg_recreate_on_change        = true
      additional_security_group_ids = [aws_security_group.worker_group_mgmt_one.id]
      labels = {
        "node_type" : "basic"
      }
    },
    {
      name                          = "worker-group-runner"
      instance_type                 = "t2.medium"
      asg_min_size                  = 1
      asg_max_size                  = 4
      asg_recreate_on_change        = true
      default_cooldown              = null
      protect_from_scale_in         = false
      additional_security_group_ids = [aws_security_group.worker_group_mgmt_two.id]
      labels = {
        "node_type" : "runner"
      }
      tags = [
        {
          "key"                 = "k8s.io/cluster-autoscaler/enabled"
          "value"               = "true"
          "propagate_at_launch" = "false"
        },
        {
          "key"                 = "k8s.io/cluster-autoscaler/${local.cluster_name}"
          "value"               = "true"
          "propagate_at_launch" = "false"
        },
        {
          "key"                 = "k8s.io/cluster-autoscaler/node-template/taint/node_type"
          "value"               = "runner:NoSchedule"
          "propagate_at_launch" = "false"
        },
        {
          "key"                 = "k8s.io/cluster-autoscaler/node-template/taint/node_type"
          "value"               = "runner:NoExecute"
          "propagate_at_launch" = "false"
        },
      ]
      kubelet_extra_args = join(" ", [
        "--node-labels=cluster_name=${local.cluster_name},node_type=runner",
        "--register-with-taints node_type=runner:NoSchedule,node_type=runner:NoExecute"
      ])
    },
  ]

  tags = {
    Environment = local.environment
    Terraform   = true
  }
}

data "aws_eks_cluster" "cluster" {
  name = module.eks.cluster_id
}

data "aws_eks_cluster_auth" "cluster" {
  name = module.eks.cluster_id
}

provider "kubernetes" {
  host                   = data.aws_eks_cluster.cluster.endpoint
  cluster_ca_certificate = base64decode(data.aws_eks_cluster.cluster.certificate_authority.0.data)
  token                  = data.aws_eks_cluster_auth.cluster.token
  load_config_file       = false
}
