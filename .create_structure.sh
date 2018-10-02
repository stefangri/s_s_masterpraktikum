read -p "projectname: " dir_name
read -p "title: " title
read -p "number: " vnr
read -p "date: " date_durchfuehrung

cp -avr .templates/structure $dir_name


echo "\newcommand{\versuch}{$title}" >> $dir_name/report/.data.tex
echo "\newcommand{\vnr}{$vnr}" >> $dir_name/report/.data.tex
echo "\newcommand{\vd}{Tag der Durchführung: $date_durchfuehrung}" >> $dir_name/report/.data.tex
echo "\newcommand{\va}{Tag der Abgabe: $date_abgabe}" >> $dir_name/report/.data.tex