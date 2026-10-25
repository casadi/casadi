function nfail = doctest_matlab(varargin)
% Run the \doctest blocks of the \lang{matlab} groups of CasADi C++ headers (MATLAB and Octave).
%
%   doctest_matlab()                      all headers under casadi/
%   doctest_matlab(file1, file2, ...)     selected headers (or directories)
%   doctest_matlab(..., '-verbose')       print every example
%
% Block syntax inside a header comment:
%   \lang{matlab}
%   \doctest
%   >> x = MX.sym('x', 2);
%   >> disp(size(x))
%        2     1
%   \enddoctest
%   \endlang
% Each '>> ' line is one statement (several separated by ';' or ',' are fine). A piece without
% semicolon shows its value like disp() does (no 'ans =' / 'x =' header, strings unquoted). The lines up to the
% next '>> ' are its expected output; blank lines and runs of blanks are ignored and a line '...'
% matches anything; numbers compare to their written precision (rule in doctest_headers.py).
% Blocks run in a fresh workspace with 'import casadi.*'.
% Returns the number of failing blocks; raises an error at the end when called with no output
% argument and something failed (so CI steps fail).

  verbose = any(strcmp(varargin, '-verbose'));
  args = varargin(~strcmp(varargin, '-verbose'));
  root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));  % docs/api/extra -> repo
  if isempty(args)
    args = {fullfile(root, 'casadi')};
  end
  files = {};
  for k = 1:numel(args)
    if isfolder(args{k})
      files = [files, doctest__find(args{k})]; %#ok<AGROW>
    else
      files{end+1} = args{k}; %#ok<AGROW>
    end
  end

  octave = exist('OCTAVE_VERSION', 'builtin') > 0;
  skip = {};
  if octave                                                   % contains(); skip list
    shims = fullfile(fileparts(mfilename('fullpath')), 'doctest_octave');
    addpath(shims);
    skiplines = regexp(fileread(fullfile(shims, 'skip.txt')), '\r?\n', 'split');
    for k = 1:numel(skiplines)
      t = regexp(skiplines{k}, '^([^#|]+)\|(.*?)\s+#', 'tokens', 'once');
      if ~isempty(t), skip(end+1, :) = {strtrim(t{1}), strtrim(t{2})}; end %#ok<AGROW>
    end
  end
  nblocks = 0; nexamples = 0; nfail = 0; nskip = 0;
  for k = 1:numel(files)
    blocks = doctest__collect(files{k});
    for b = 1:numel(blocks)
      ex = blocks{b}.examples;
      [~, hdr, ext] = fileparts(files{k});
      listed = ~isempty(skip) && any(strcmp(skip(:, 1), [hdr ext]) & ismember(skip(:, 2), strtrim({ex.src})));
      if octave && (listed || any(~cellfun(@isempty, regexp({ex.src}, '<\s*casadi\.Callback', 'once'))))
        nskip = nskip + 1;                                    % Octave: MATLAB-only block
        continue
      end
      nblocks = nblocks + 1;
      nexamples = nexamples + numel(ex);
      got = doctest__run({ex.src});
      ok = true;
      for e = 1:numel(ex)
        pass = doctest__match(ex(e).want, got{e});
        if verbose || ~pass
          fprintf('%s (line %d)\n  >> %s\n', blocks{b}.where, ex(e).line, ex(e).src);
        end
        if ~pass
          ok = false;
          fprintf('Expected:\n%s\nGot:\n%s\n', doctest__indent(ex(e).want), doctest__indent(got{e}));
        end
      end
      if ~ok
        nfail = nfail + 1;
      end
    end
  end
  fprintf('matlab: %d blocks, %d examples, %d failing blocks', nblocks, nexamples, nfail);
  if nskip > 0
    fprintf(', %d skipped in Octave (docs/api/extra/doctest_octave/skip.txt, Callback subclasses)', nskip);
  end
  fprintf('\n');
  if nargout == 0
    if nfail > 0
      error('doctest_matlab:failed', '%d failing doctest block(s)', nfail);
    end
    clear nfail
  end
end

function files = doctest__find(dirname)
  files = {};
  d = dir(dirname);
  for k = 1:numel(d)
    if d(k).name(1) == '.' || (d(k).isdir && strcmp(d(k).name, 'include')), continue; end  % installed header copies
    p = fullfile(dirname, d(k).name);
    if d(k).isdir
      files = [files, doctest__find(p)]; %#ok<AGROW>
    elseif numel(p) > 4 && strcmp(p(end-3:end), '.hpp')
      files{end+1} = p; %#ok<AGROW>
    end
  end
end

function blocks = doctest__collect(fname)
  % Return a cell array of structs with fields where and examples(src, want, line).
  txt = fileread(fname);
  lines = regexp(txt, '\r?\n', 'split');
  blocks = {};
  inblock = false;
  lang = '';
  for n = 1:numel(lines)
    s = regexprep(lines{n}, '^\s*\*?\s?', '');
    t = strtrim(regexprep(s, '\*/\s*$', ''));  % '\enddoctest */' closes the comment too
    l = regexp(t, '^\\lang\{(\w+)\}$', 'tokens', 'once');
    if ~inblock && ~isempty(l)
      lang = l{1};
    elseif ~inblock && strcmp(t, '\endlang')
      lang = '';
    elseif ~inblock && strcmp(t, '\doctest') && strcmp(lang, 'matlab')
      inblock = true; start = n; body = {}; bodyline = [];
    elseif inblock && strcmp(t, '\enddoctest')
      inblock = false;
      blocks{end+1} = struct('where', sprintf('%s:%d', fname, start), ...
        'examples', doctest__examples(body, bodyline)); %#ok<AGROW>
    elseif inblock
      body{end+1} = s; %#ok<AGROW>
      bodyline(end+1) = n; %#ok<AGROW>
    end
  end
  if inblock
    error('doctest_matlab:syntax', '%s:%d: unterminated \\doctest', fname, start);
  end
end

function ex = doctest__examples(body, bodyline)
  ex = struct('src', {}, 'want', {}, 'line', {});
  for n = 1:numel(body)
    if numel(body{n}) >= 3 && strcmp(body{n}(1:3), '>> ')
      ex(end+1) = struct('src', body{n}(4:end), 'want', '', 'line', bodyline(n)); %#ok<AGROW>
    elseif ~isempty(ex)
      if isempty(ex(end).want)
        ex(end).want = body{n};
      else
        ex(end).want = sprintf('%s\n%s', ex(end).want, body{n});
      end
    end
  end
end

function doctest__outs = doctest__run(doctest__stmts)
  % Runs in its own workspace: variables of one block never leak into the next.
  import casadi.*
  doctest__outs = cell(size(doctest__stmts));
  for doctest__k = 1:numel(doctest__stmts)
    doctest__s = doctest__stmts{doctest__k};
    doctest__o = '';
    try
      if ~isempty(regexp(doctest__s, '^\s*(try|if|for|parfor|while|switch|function)(?=\W|$)', 'once'))
        doctest__o = evalc(doctest__s);                       % control flow: as written
      else
        doctest__p = doctest__pieces(doctest__s);
        for doctest__j = 1:size(doctest__p, 1)
          doctest__e = doctest__p{doctest__j, 1};
          if doctest__p{doctest__j, 2}                        % ';': no display
            doctest__o = [doctest__o evalc([doctest__e ';'])]; %#ok<AGROW>
          elseif doctest__is_assignment(doctest__e)           % 'x = ...': value without header
            doctest__o = [doctest__o evalc([doctest__e '; disp(' doctest__lhs(doctest__e) ')'])]; %#ok<AGROW>
          elseif ~isempty(regexp(doctest__e, '^\s*((disp|fprintf|display|print)\s*\(|(cellfun|arrayfun)\s*\(\s*@\(?\w*\)?\s*(disp|fprintf))', 'once'))
            doctest__o = [doctest__o evalc(doctest__e)]; %#ok<AGROW>
          else                                                % bare expression: shown like disp
            doctest__rerun = false;
            doctest__x = evalc(['try, disp(' doctest__e '), catch doctest__ME, ' ...
              'doctest__rerun = doctest__void(doctest__ME); end']);
            if doctest__rerun                                 % rejected before running: run as is
              doctest__x = evalc(doctest__e);
            end
            doctest__o = [doctest__o doctest__x]; %#ok<AGROW>
          end
        end
      end
      doctest__outs{doctest__k} = doctest__o;
    catch doctest__err
      doctest__outs{doctest__k} = [doctest__o 'ERROR: ' doctest__err.message];
    end
  end
end

function rerun = doctest__void(ME)
  % disp(f(...)) of a call without return value. 'Too many output arguments' is raised before
  % the call runs (run it again without disp); 'not assigned' after it ran (keep its output).
  if ~isempty(regexp(ME.message, 'Too many output arguments', 'once'))
    rerun = true;
  elseif ~isempty(regexp(ME.message, '(output arguments? not assigned|value expected|undefined in return list)', 'once'))
    rerun = false;
  else
    rethrow(ME);
  end
end

function p = doctest__pieces(s)
  % Split a statement line on top-level ',' and ';' (outside brackets and string literals).
  % Rows: {text, suppressed}.
  p = cell(0, 2); depth = 0; inq = false; cur = ''; prev = ' '; i = 0;
  while i < numel(s)
    i = i + 1; c = s(i);
    if inq
      cur = [cur c]; %#ok<AGROW>
      if c == '''' && i < numel(s) && s(i + 1) == ''''       % '' inside a string literal
        cur = [cur ''''];  i = i + 1; %#ok<AGROW>
      elseif c == ''''
        inq = false;
      end
    elseif c == '''' && isempty(regexp(prev, '[\w\)\]\}\.'']', 'once'))
      inq = true; cur = [cur c]; %#ok<AGROW>
    elseif any(c == '([{')
      depth = depth + 1; cur = [cur c]; %#ok<AGROW>
    elseif any(c == ')]}')
      depth = depth - 1; cur = [cur c]; %#ok<AGROW>
    elseif (c == ',' || c == ';') && depth == 0
      if ~isempty(strtrim(cur)); p(end+1, :) = {strtrim(cur), c == ';'}; end %#ok<AGROW>
      cur = '';
    else
      cur = [cur c]; %#ok<AGROW>
    end
    if ~isspace(c); prev = c; end
  end
  if ~isempty(strtrim(cur)); p(end+1, :) = {strtrim(cur), false}; end
end

function tf = doctest__is_assignment(e)
  tf = ~isempty(regexp(e, '^\s*(\[[^\]]*\]|[A-Za-z]\w*(\.\w+|\([^)]*\)|\{[^}]*\})*)\s*=(?!=)', 'once'));
end

function v = doctest__lhs(e)
  % Variables to show after an unsuppressed assignment: 'x' or each name of '[a, b]'.
  t = regexp(e, '^\s*(\[[^\]]*\]|[A-Za-z]\w*)', 'tokens', 'once');
  v = strtrim(t{1});
  if v(1) == '['
    v = strjoin(regexp(v(2:end-1), '[A-Za-z]\w*', 'match'), '), disp(');
  end
end

function ok = doctest__match(want, got)
  w = doctest__norm(want);
  g = doctest__as_expected(w, doctest__norm(got));
  if strcmp(w, g)
    ok = true;
    return
  end
  if isempty(strfind(w, '...'))
    ok = false;
    return
  end
  parts = strsplit(w, '...');
  for k = 1:numel(parts)
    parts{k} = regexptranslate('escape', parts{k});
  end
  ok = ~isempty(regexp(g, ['^' strjoin(parts, '.*') '$'], 'once'));
end

function g = doctest__as_expected(w, g)
  % Numbers as in doctest_headers.py: half a unit in the last written digit, at least 1e-9 relative.
  pat = '(?<![\w.@])[-+]?(\d+\.\d*|\.\d+|\d+)([eE][-+]?\d+)?(?![\w.])';
  wm = regexp(w, pat, 'match');
  wm = wm(cellfun(@isempty, regexp(wm, '^[-+]?0\d', 'once')));
  [gs, ge, gm] = regexp(g, pat, 'start', 'end', 'match');
  keep = cellfun(@isempty, regexp(gm, '^[-+]?0\d', 'once'));
  gs = gs(keep); ge = ge(keep); gm = gm(keep);
  if numel(wm) ~= numel(gm)
    return
  end
  for k = numel(gm):-1:1
    if abs(str2double(wm{k}) - str2double(gm{k})) <= doctest__tolerance(wm{k}) * (1 + 1e-9)
      g = [g(1:gs(k)-1) wm{k} g(ge(k)+1:end)];
    end
  end
end

function t = doctest__tolerance(lit)
  t = 0;
  if any(lit == '.') || any(lit == 'e') || any(lit == 'E')
    d = regexp(lit, '\.(\d*)', 'tokens', 'once');
    x = regexp(lit, '[eE]([-+]?\d+)', 'tokens', 'once');
    nd = 0; ne = 0;
    if ~isempty(d), nd = numel(d{1}); end
    if ~isempty(x), ne = str2double(x{1}); end
    t = 0.5 * 10^(ne - nd);
  end
  x = abs(str2double(lit));
  if x == 0, x = 1; end
  t = max(t, 1e-9 * x);
end

function s = doctest__norm(s)
  lines = regexp(s, '\r?\n', 'split');
  keep = {};
  for k = 1:numel(lines)
    l = strtrim(regexprep(lines{k}, '\s+', ' '));
    if ~isempty(l)
      keep{end+1} = l; %#ok<AGROW>
    end
  end
  s = strjoin(keep, char(10));
end

function s = doctest__indent(s)
  s = ['    ' strrep(strtrim(s), char(10), [char(10) '    '])];
end
