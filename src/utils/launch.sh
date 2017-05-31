#!/bin/bash
for i in {1..10}
do
		./fastnn allfeatures$i $i> test$i.txt &
done

