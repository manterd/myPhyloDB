#!/bin/sh

echo "Installing myPhyloDB vers. 1.1\n"
echo "Checking if myPhyloDB exists...\n"

if [ -d "$HOME/myPhyloDB" ]
    then
	echo "A previous installation of myPhyloDB was detected..."
        echo -n "Would you like to keep your current database files [y/n]?"
        read response
    else
	echo "No previous installations of myPhyloDB were detected..."
        response=n
fi

if [ $response = y ]
    then
	mkdir $HOME/myPhyloDB_temp
	cp $HOME/myPhyloDB/dbMicrobe $HOME/myPhyloDB_temp/dbMicrobe
	cp -r $HOME/myPhyloDB/uploads $HOME/myPhyloDB_temp/uploads
	tar -zxf myPhyloDB.tar.gz -C $HOME
	rm $HOME/myPhyloDB/dbMicrobe
	rm -rf $HOME/myPhyloDB/uploads
	cp $HOME/myPhyloDB_temp/dbMicrobe $HOME/myPhyloDB/dbMicrobe
	cp -r $HOME/myPhyloDB_temp/uploads $HOME/myPhyloDB/uploads
	rm -rf $HOME/myPhyloDB_temp
    else
	if [ -d "$HOME/myPhyloDB" ]
	    then
		rm -rf $HOME/myPhyloDB
	fi
	tar -zxf myPhyloDB.tar.gz -C $HOME
fi

mkdir -p $HOME/.icons
cp $HOME/myPhyloDB/media/images/myPhyloDB_Logo.png $HOME/.icons/myPhyloDB_Logo.png
cp myPhyloDB.desktop $HOME/Desktop/myPhyloDB.desktop
chmod +x $HOME/Desktop/myPhyloDB.desktop

echo ""
echo "myPhyloDB v.1.1 installation is finished!"
