tar czvf ath.tgz athena-graphics-3d
tar czvf par.tgz par
docker build --tag mhd .
docker run -it mhd /bin/bash
