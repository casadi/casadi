function [out]=casadi_init_struct()
    try
        out = struct();
        out.defaultCompiler = '';
        out.defaultLinker = '';
        out.defaultCompilerSetup = '';
        out.defaultLinkerSetup = '';
        out.defaultCompilerOutputFlag = '';
        out.defaultLinkerOutputFlag = '';
    catch ME
        out = struct();
        out.error = getReport(ME);
    end    
end
