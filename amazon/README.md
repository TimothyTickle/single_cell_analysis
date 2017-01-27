# Amazon setup

## pull single cell demo code from github

    git clone https://github.com/TimothyTickle/single_cell_analysis.git

## pull docker image from docker hub

    docker pull trinityctat/cegs_demo


# Launch single container for testing

    docker run --rm -p 8787:8787 -v `pwd`/single_cell_analysis/src:/home/training/src:ro --name my_rstudio_docker trinityctat/cegs_demo




## Misc notes on Amazon EC2 setup

----
#!/bin/bash
curl https://get.docker.com | sh
sudo usermod -a -G docker ubuntu
sudo service docker start
----

