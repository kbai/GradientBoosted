#!/bin/bash
for i in {1..10}
do
		./fastnn test$i > test$i.txt &
done

