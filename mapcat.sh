#!/bin/bash

if [ $2 -lt 10 ] ; then
	if [ $3 -lt 9 ] ; then
		for (( i=$2 ; i <= $3 ; i++ ))
		do
			cat $1*_000${i}_phi.map > $1_000${i}_phi.map
		done
	else
		for (( i=$2 ; i <= 9 ; i++ ))
		do
			cat $1*_000${i}_phi.map > $1_000${i}_phi.map
		done
		if [ $3 -lt 99 ] ; then
			for (( i=10 ; i <= $3 ; i++ ))
			do
				cat $1*_00${i}_phi.map > $1_00${i}_phi.map
			done
		else
			for (( i=10 ; i <= 99 ; i++ ))
			do
				cat $1*_00${i}_phi.map > $1_00${i}_phi.map
			done
			if [ $3 -lt 999 ] ; then
				for (( i=100 ; i <= $3 ; i++ ))
				do
					cat $1*_0${i}_phi.map > $1_0${i}_phi.map
				done
			else
				for (( i=100 ; i <= 999 ; i++ ))
				do
					cat $1*_0${i}_phi.map > $1_0${i}_phi.map
				done
				for (( i=1000 ; i <= $3 ; i++ ))
				do
					cat $1*_${i}_phi.map > $1_${i}_phi.map
				done
			fi
		fi
	fi
else
	if [ $2 -lt 100 ] ; then
		if [ $3 -lt 99 ] ; then
			for (( i=$2 ; i <= $3 ; i++ ))
			do
				cat $1*_00${i}_phi.map > $1_00${i}_phi.map
			done
		else
			for (( i=$2 ; i <= 99 ; i++ ))
			do
				cat $1*_00${i}_phi.map > $1_00${i}_phi.map
			done
			if [ $3 -lt 999 ] ; then
				for (( i=100 ; i <= $3 ; i++ ))
				do
					cat $1*_0${i}_phi.map > $1_0${i}_phi.map
				done
			else
				for (( i=100 ; i <= 999 ; i++ ))
				do
					cat $1*_0${i}_phi.map > $1_0${i}_phi.map
				done
				for (( i=1000 ; i <= $3 ; i++ ))
				do
					cat $1*_${i}_phi.map > $1_${i}_phi.map
				done
			fi
		fi
	else
		if [ $3 -lt 999 ] ; then
			for (( i=$2 ; i <= $3 ; i++ ))
			do
				cat $1*_0${i}_phi.map > $1_0${i}_phi.map
			done
		else
			for (( i=$2 ; i <= 999 ; i++ ))
			do
				cat $1*_0${i}_phi.map > $1_0${i}_phi.map
			done
			for (( i=1000 ; i <= $3 ; i++ ))
			do
				cat $1*_${i}_phi.map > $1_${i}_phi.map
			done
		fi
	fi
fi
exit 0
