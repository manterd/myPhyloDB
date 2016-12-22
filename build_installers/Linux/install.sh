#!/bin/sh

# This file will not install myPhyloDB.
# To install, please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472


echo "Installing myPhyloDB v.1.2.0\n"
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

mkdir $HOME/myPhyloDB_temp
mkdir $HOME/myPhyloDB_temp/uploads

if [ $response = y ]
    then
        mv $HOME/myPhyloDB/db.Microbe $HOME/myPhyloDB_temp/db.Microbe
        mv $HOME/myPhyloDB/db.PICRUSt $HOME/myPhyloDB_temp/db.PICRUSt
        mv $HOME/myPhyloDB/uploads/* $HOME/myPhyloDB_temp/uploads/
        rm -rf $HOME/myPhyloDB
        tar -zxf myPhyloDB.tar.gz -C $HOME
        mv -f $HOME/myPhyloDB_temp/db.Microbe $HOME/myPhyloDB/db.Microbe
        mv -f $HOME/myPhyloDB_temp/db.PICRUSt $HOME/myPhyloDB/db.PICRUSt
        rm -rf $HOME/myPhyloDB/uploads
        mkdir $HOME/myPhyloDB/uploads
        mv -f $HOME/myPhyloDB_temp/uploads/* $HOME/myPhyloDB/uploads/
    else
        if [ -d "$HOME/myPhyloDB" ]
            then
                mv $HOME/myPhyloDB/uploads/* $HOME/myPhyloDB_temp/uploads/
                rm -rf $HOME/myPhyloDB
        fi
        tar -zxf myPhyloDB.tar.gz -C $HOME
        if [ -d "$HOME/myPhyloDB_temp/uploads" ]
            then
                rm -rf $HOME/myPhyloDB/uploads
                mkdir $HOME/myPhyloDB/uploads
                mv -f $HOME/myPhyloDB_temp/uploads/* $HOME/myPhyloDB/uploads/
        fi
fi

rm -rf $HOME/myPhyloDB_temp

mkdir -p $HOME/.icons
cp $HOME/myPhyloDB/media/images/myPhyloDB_Logo.png $HOME/.icons/myPhyloDB_Logo.png
cp myPhyloDB.desktop $HOME/Desktop/myPhyloDB.desktop
chmod +x $HOME/Desktop/myPhyloDB.desktop

# fix path in R file
sed -i 's/PycharmProjects\///g' $HOME/myPhyloDB/R/R-Linux/bin/R
sed -i 's/PycharmProjects\///g' $HOME/myPhyloDB/R/R-Linux/lib/R/bin/R

echo ""
echo "myPhyloDB v.1.2.0 installation is finished!"
