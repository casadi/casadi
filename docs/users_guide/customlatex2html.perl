package main;

sub do_cmd_verbatiminputeval {
  local($outer) = @_;
  $outer=~ s/\\arabic(<<\d+>>\w+<<\d+>>)/(&read_counter_value($1))[1]/e;
  $outer=~ s/\\alph(<<\d+>>\w+<<\d+>>)/&falph((&read_counter_value($1))[1])/e;
  $outer=~ s/\\verbatiminputeval/\\verbatiminput/;
  &do_cmd_verbatiminput($outer);
}

&process_commands_wrap_deferred (<<_RAW_ARG_DEFERRED_CMDS_);
verbatiminputeval # {}
_RAW_ARG_DEFERRED_CMDS_

1;
