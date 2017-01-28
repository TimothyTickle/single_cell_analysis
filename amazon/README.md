# Amazon setup

## Launch an instance:

*  From <https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#Instances:sort=instanceId>, launch an EC2 instance.
*  Step 1: Select 64-bit Ubuntu server 16.04
*  Step 2: select instance type 't2.large' (2 vcpus, 8G RAM), click button 'Next: Configure Instance Details'
*  Step 3: Click 'Advanced Details' at bottom. In the text box shown, insert the following text:

    #!/bin/bash
    curl https://get.docker.com | sh
    sudo usermod -a -G docker ubuntu
    sudo service docker start

then, click 'Next: Add Storage'

* Step 4: 'Storage' set size to 25

* Step 5: 'Add tags' - skip this step.  Click button 'Next: Configure Security Group'

* Step 6: 'Configure Security Group': Click 'Add Rule', and add:

    Custom TCP Rule      TCP      8787      Anywhere

* You're done.  Click 'Review and Launch', followed by Launch.

It'll ask you to select or create a new key pair.  Assuming you don't already have one, select 'Create a new key pair', and then download that keypair to your local computer.  Click 'Launch Instances'.

Click 'View Instances' at the bottom right of the page.   This will bring you to the EC2 dashboard.  You should see your new instance listed in the dashboard, and the column 'IPv4 Public IP' will contain the IP address that you'll use to connect to your instance using ssh to complete the configuration.

## Configure your EC2 Instance for the Single Cell Workshop

SSH to your instance using your key pair .pem file:

    ssh -i your_ssh_key.pem ubuntu@${IPv4_Public_IP}

using the IP address from the corresponding column in the EC2 dashboard as described above.

Complete the following steps below to pull down the code and data from github, pull down the Docker container, and launch the docker container.

### pull single cell demo code from github

    git clone https://github.com/TimothyTickle/single_cell_analysis.git

### pull docker image from docker hub

    docker pull trinityctat/cegs_demo


### Launch single container for testing

    docker run --rm -p 8787:8787 -v `pwd`/single_cell_analysis/src:/home/training/src:ro --name my_rstudio_docker -d trinityctat/cegs_demo



## Run RStudio

From within your web browser, go to the link http://${IPv4_Public_IP}:8787 using your IP address, and the ':8787' indicates that it will connect to the port that RStudio listens to on your instance.

Assuming RStudio runs within your browser, it may prompt you for user/password, in which case you can use 'training/training'.

Now you're ready to run through the workshop exercises.  Good luck!

