echo "Graph path: "
read graphTxt

for k in `seq 2 3`
do
	./cn-order $graphTxt $k >> "web_google_cn_time.txt"
done
