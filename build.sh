#!/bin/sh
gcc -g main.c linalg.c -I./ -fsanitize=address,undefined -Wall -Werror -o linalg
