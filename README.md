# Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This repository contains all the necessary code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data, available in 

* C Soneson & MD Robinson: [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612). Nature Methods 15:255-261 (2018).


## Docker container

Build using

```bash
docker build -t conquer .
```

from the *conquer_comparison* root directory.

Run the container with 

```bash
docker run --rm \
    -e PASSWORD=conquer \
    -p 8787:8787 \
    -v /Users/milan/Projects/conquer_test:/home/conquer \
    conquer
```

and access RStudio instance at `localhost:8787`. `--rm` removes the container automatically once it exits.

To access the running container via command line, run

```
docker exec -it <container ID> /bin/bash
```

__after__ the container has been launched. This should open a bash shell inside the container.
