W=$1
tmux new -d -s $W nice ./submit.sh $W
