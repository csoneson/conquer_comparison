# Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This repository contains all the necessary code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data, available in 

> C Soneson & MD Robinson: [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612). Nature Methods 15:255-261 (2018).


## Docker container

Run the container with 

```bash
docker run --rm -d \
    -e PASSWORD=conquer \
    -v /Users/milan/Projects/conquer_test:/home/conquer \
    conquer
```

To access the running container via command line, run

```
docker exec -it <container ID> /bin/bash
```

__after__ the container has been launched. This should open a bash shell inside the container. The container ID can be found by checking the active containers with `docker ps -a`.

__Note__: some of the steps in the pipeline require quite some memory, so make sure to assign __at least 4GB__ of memory to your Docker instance (default is 2GB). Instructions for this can be found [here](https://stackoverflow.com/a/44533437/11801854).
