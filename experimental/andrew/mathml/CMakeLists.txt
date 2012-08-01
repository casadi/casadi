add_custom_command(OUTPUT example1.xml COMMAND matlab -nodesktop -nojvm -nosplash -r 'run examples; exit' DEPENDS examples.m)
