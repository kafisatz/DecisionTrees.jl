cmd = `$(s.graphvizexecutable) "-V"`
gres = run(cmd)

(path,io) = mktemp()
run(gcmd ,stdout=io)

run(pipeline(gcmd,stdout=io));

tmptxt = mktemp()[1]
tmperr = mktemp()[1]
pip = pipeline(cmd, stdout=tmptxt, stderr=tmperr)
rc = run(pip)
        