% This is a _very_ primitive cross-platform patch utility

% === begin Setup =====
orig_name = "octrunORIG.swg";  % The input file you want to see patched
result_name = "octrun.swg";    % The file of the patched result. Should not be equal to orig_name. 

% Serial stack of patches to be applied in (match_string_i, filename_i) style.
% If match_string_i is found in the orig_name file, the contents of filename_i will be inserted before this line
% After this, i is increased by one.
patches = {
{"virtual bool is_string() const", "octave_swig_type.snippet"}; 
{"virtual bool is_string() const", "octave_swig_ref.snippet"};
};
% === End Setup ====




disp(["This is octave patcher, working on " orig_name " => " result_name])

orig    = fopen(orig_name,  "r");
if orig==-1
  error(["file '" orig_name "' not found in '" pwd "'"]);
end

patched = fopen(result_name,"w");
try

i = 1;
orig_line = fgetl(orig);
while ischar(orig_line)
    % Check if there are patches on the serial stack
    if i <= length(patches)
      if index(orig_line,patches{i}{1})
        disp(["found snippet '" patches{i}{1} "' => '" patches{i}{2} "'"])
        snip = fopen(patches{i}{2},"r");
        if snip==-1
          error(["file '" patches{i}{2} "' not found in '" pwd "'"])
        end
        patch_line = fgetl(snip);
        while ischar(patch_line)
          fdisp(patched,patch_line);
          patch_line = fgetl(snip);
        end
        fclose(snip);
        i = i + 1;
      end
    end
    % Copy line from original file to result file
    fdisp(patched,orig_line);
    orig_line  = fgetl(orig);
end

fclose(orig);
fclose(patched);
disp(["Octave patcher completed succesfully."])

catch
  try
     fclose(patched);
  end_try_catch      
  unlink(result_name);  % Don't leave behind an interrupted result file
  disp(lasterr)
  disp("Octave patcher failed. Falling back to unpatched SWIG");
end_try_catch
