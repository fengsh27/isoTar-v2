#!/bin/bash
# Kill all the damn Zombie processes!

kill $(ps -A -ostat,ppid | awk '/[zZ]/{print $2}')
