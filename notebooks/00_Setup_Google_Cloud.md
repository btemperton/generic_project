# Using GoogleCloud for Snakemake

Using GoogleCloud, it is possible to run your jobs on the Google infrastructure for a fee. It does this via Kubernetes.

## Before you begin
Kubernetes uses the `mb_mem` and `threads` flags to determine which nodes to send jobs to in your Kubernetes cluster, so make sure they are all set in your rules.

## Set up Bucket for Storage

```
PROJECT_NAME=wec-test
ZONE=europe-west2
STORAGE_NAME=tempertonlab-wec-store
gsutil mb -c standard -l $ZONE -p $PROJECT_NAME gs://$STORAGE_NAME/

```

## Set up Kubernetes Cluster

```
CLUSTER_NAME=wec-cluster
ZONE=europe-west2-a

gcloud container clusters create $CLUSTER_NAME \
--zone $ZONE \
--scopes storage-rw \
--machine-type=n1-highmem-96 \
--enable-autoscaling \
--min-nodes=0 \
--max-nodes=16 \
--num-nodes=1

gcloud container clusters get-credentials $CLUSTER_NAME --zone=$ZONE
```

##Running Snakemake with Kubernetes

```
snakemake \
-p \
--kubernetes \
--use-conda \
-j 12 \
--default-remote-provider GS \
--default-remote-prefix $STORAGE_NAME
```
