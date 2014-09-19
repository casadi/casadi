while (<>) {
  $description = $1 if /\* (.*) \*/;
  $type = "OT_BOOL" if /^\s+bool/;
  $type = "OT_REAL" if /^\s+double/;
  $type = "OT_INTEGER" if /^\s+int/;

  if (/^\s+(\w+)\s+(\w+)/ and $2 ne"Ares") {
    print qq|addOption("$2",$type,worhp_p.$2,"$description");\n|;
    #print qq|if (hasSetOption("$2")) worhp_p.$2 = getOption("$2");\n|;
  }

}
