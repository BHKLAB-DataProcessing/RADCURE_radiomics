
Set env variables for dockerfile:
``` bash
export DOCKER_REPO="jjjermiah"
export DOCKER_IMAGE_NAME="radcure_radiomics"
export DOCKER_TAG="0.1"
```


Commands to build the docker container:
``` bash

docker build -t $DOCKER_REPO/$DOCKER_IMAGE_NAME:$DOCKER_TAG -f build/R_steps.Dockerfile .
```

Command to push to dockerhub:
``` bash
docker push $DOCKER_REPO/$DOCKER_IMAGE_NAME:$DOCKER_TAG
```