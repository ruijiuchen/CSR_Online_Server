while [ 1 ]
do
    online_server=`pgrep CSR_Online_Serv`
    echo $online_server
    if [ -n "$online_server" ]
    then
        echo "online_server pid: $online_server"
	kill -10 $online_server
	sleep 5
    fi
    sleep 1
done
