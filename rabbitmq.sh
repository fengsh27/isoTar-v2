#!/bin/bash
/bin/sleep 4
/usr/sbin/rabbitmqctl add_user isorabmq isorabmq
/usr/sbin/rabbitmqctl add_vhost isorabmq_vhost
/usr/sbin/rabbitmqctl set_user_tags isorabmq isorabmq_tag
/usr/sbin/rabbitmqctl set_permissions -p isorabmq_vhost isorabmq ".*" ".*" ".*"
service rabbitmq-server restart
