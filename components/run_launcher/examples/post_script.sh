# Post-script example
# NOTE: You can access $run variables with $run.sub_value.another_sub_value.etc
# this shouldn't have issues as mongoDB doesn't allow .'s in key names so we can
# parse on that value
echo "start post_script $run.name $run.type";
echo "end post_script";
