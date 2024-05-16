# Setup

First, open a the terminal (mac/linux) or the command prompt (windows), and navigate to the folder you would like to add the shoreline change model directory to. Then clone the repo from github:

```sh
git clone https://github.com/Shannon-Bengtson/PacificCoastalHazards.git
```

Setup the environment using conda.

```sh
conda env create -f environment.yml -n myenv
```

And then activte the environment

```sh
conda activate myenv
```

If running locally on your desktop, we can simply open jupyter lab and run the notebooks from there.

```sh
jupyter lab
```

If you are running on an external server via ssh, we need to create tunnel through which to open the notebooks. On the external server (where the environment and github repo are located), run:

```sh
jupyter lab --no-browser --port=8989
```

Then, open a new terminal/command prompt on your local computer and run the following, replacing user@hostname with your username and the server IP address or domain:

```sh
ssh -N -f -L localhost:8989:localhost:8989 <user>@hostname
```

Then, in a browser naviate to localhost:8989. If a token is required, copy and paste the token key from the output in terminal/command prompt window for the external server.
