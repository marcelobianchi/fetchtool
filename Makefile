all:
	echo "Choose one:"
	echo ""
	echo "make pipinstall"
	echo "make pipremove"
	echo ""

pipinstall:
	pip3 install -e --no-deps ./

pipremove:
	pip3 uninstall --yes fetchtool

