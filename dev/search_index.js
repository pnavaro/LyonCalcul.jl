var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#LyonCalcul.Advection",
    "page": "Documentation",
    "title": "LyonCalcul.Advection",
    "category": "type",
    "text": "Advection(n, p, delta)\n\nn :: Number of points.\np :: Spline degree.\ndelta :: space size step.\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.Advection-Tuple{Array{Float64,1},Float64}",
    "page": "Documentation",
    "title": "LyonCalcul.Advection",
    "category": "method",
    "text": "advection! = Advection( n, p, delta)\nadvection!( f, alpha )\n\nCreate a function to compute the interpolating spline \nof degree p of odd\ndegree of a 1D function f on a periodic uniform mesh, at\nall points x after a displacement alpha. \nInput f type is Vector{Float64} and is updated inplace.\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.bspline-Tuple{Int64,Int64,Float64}",
    "page": "Documentation",
    "title": "LyonCalcul.bspline",
    "category": "method",
    "text": "bspline(p, j, x)\n\nReturn the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j. Implemented recursively using the De Boor\'s Algorithm\n\nB_i0(x) = left\nbeginmatrix\n1  mathrmif  quad t_i  x  t_i+1 \n0  mathrmotherwise \nendmatrix\nright\n\nB_ip(x) = fracx - t_it_i+p - t_i B_ip-1(x) \n+ fract_i+p+1 - xt_i+p+1 - t_i+1 B_i+1p-1(x)\n\n\n\n\n\n"
},

{
    "location": "#LyonCalcul.jl-1",
    "page": "Documentation",
    "title": "LyonCalcul.jl",
    "category": "section",
    "text": "Documentation for LyonCalcul.jlModules = [LyonCalcul]\nOrder   = [:type, :function]"
},

{
    "location": "slides/#",
    "page": "Talk",
    "title": "Talk",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\""
},

{
    "location": "slides/#Who-am-I-?-1",
    "page": "Talk",
    "title": "Who am I ?",
    "category": "section",
    "text": "My name is Pierre Navaro\nPh.D in Computational Aeroacoustics, 1998-2002 (Université du Havre) (Fortran 77+PVM)\nScientific Software Engineer in Strasbourg (2003-2015) (Fortran 90-2003 + OpenMP-MPI)\nMoved to Rennes in 2015 (Numpy + Cython, R + Rcpp)\nJulia user since July 2018 (Julia v1.0)"
},

{
    "location": "slides/#bspline-function-1",
    "page": "Talk",
    "title": "bspline function",
    "category": "section",
    "text": "\"\"\"\n    bspline(p, j, x)\n\nReturn the value at x in [0,1[ of the B-spline with integer nodes\nof degree p with support starting at j.\nImplemented recursively using the\n[De Boor\'s Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)\nmath B{i,0}(x) := \\left\\{ \\begin{matrix} 1 & \\mathrm{if}  \\quad ti ≤ x < t_{i+1} \\\\\n0 & \\mathrm{otherwise} \\end{matrix} \\right.math B{i,p}(x) := \\frac{x - ti}{t{i+p} - ti} B_{i,p-1}(x)\\frac{t{i+p+1} - x}{t{i+p+1} - t{i+1}} B{i+1,p-1}(x).\"\"\"\nfunction bspline(p::Int, j::Int, x::Float64)\n   if p == 0\n       if j == 0\n           return 1.0\n       else\n           return 0.0\n       end\n   else\n       w = (x - j) / p\n       w1 = (x - j - 1) / p\n   end\n   return (w * bspline(p - 1, j, x)\n           + (1 - w1) * bspline(p - 1, j + 1, x))\nend?bspline"
},

{
    "location": "slides/#Advection-callable-struct-1",
    "page": "Talk",
    "title": "Advection callable struct",
    "category": "section",
    "text": "\"\"\"\n    Advection(n, p, delta)\n\n    n :: Number of points.\n    p :: Spline degree.\n    delta :: space size step.\n\n\"\"\"\nmutable struct Advection\n    n         :: Int64\n    p         :: Int64\n    delta     :: Float64\n    modes     :: Vector{Float64}\n    eig_bspl  :: Vector{ComplexF64}\n    eig_alpha :: Vector{ComplexF64}\n\n    function Advection( n, p, delta )\n\n        eig_bspl  = zeros(ComplexF64,n)\n        eig_alpha = similar(eig_bspl)\n        modes = 2π * (0:n-1) / n\n\n        eig_bspl .= bspline(p, -div(p+1,2), 0.0)\n        for j in 1:div(p+1,2)-1\n            eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)\n             * 2 * cos.(j * modes))\n        end\n\n        new(n, p, delta, modes, eig_bspl, eig_alpha)\n    end\nendusing FFTW\n\n\"\"\"\n    advection! = Advection( n, p, delta)\n    advection!( f, alpha )\n\n    Create a function to compute the interpolating spline\n    of degree p of odd degree of a 1D function f on a periodic\n    uniform mesh, at all points x after a displacement alpha.\n    Input f type is Vector{Float64} and is updated inplace.\n\n\"\"\"\nfunction (adv :: Advection)(f     :: Vector{Float64},\n                            alpha :: Float64)\n\n\n   ishift = floor(- alpha / adv.delta)\n   beta = - ishift - alpha / adv.delta\n   fill!(adv.eig_alpha, 0.0)\n   for j in -div(adv.p-1,2):div(adv.p+1,2)\n      adv.eig_alpha .+= (bspline(adv.p, j-div(adv.p+1,2), beta)\n         .* exp.((ishift + j) * 1im .* adv.modes))\n   end\n\n   f .= real(ifft(fft(f) .* adv.eig_alpha ./ adv.eig_bspl))\n\nendusing Plots\nxmin, xmax, nx = 0, 4π, 64  ## set-up the 1d mesh\ndx = (xmax - xmin) / nx\nx = range(xmin, stop=xmax, length=nx+1)[1:end-1] |> collect\n\nf = sin.(x) ## initialize the function f(x)\nadvection! = Advection(nx, 5, dx) ## Create the advection function\nadvection!(f, 1.0) ## advect the function f with a displacement α = 1.0\nf_ref = sin.(x.-1) ## compute the analytical solution\nplot(x, f_ref, label=:reference)\nscatter!(x, f, label=:computed)"
},

{
    "location": "slides/#Why-make-a-package-?-1",
    "page": "Talk",
    "title": "Why make a package ?",
    "category": "section",
    "text": "Share your code and maintain it.\nBetter source files organization.\nImprove your programming practices (tests and documentation).\nContinuous integration (Github+Travis-CI).\nDocumentation is hosted and generated after every changes."
},

{
    "location": "slides/#Configure-git-1",
    "page": "Talk",
    "title": "Configure git",
    "category": "section",
    "text": "git config --global user.name \"Pierre Navaro\"\ngit config --global user.email \"pierre.navaro@math.cnrs.fr\"\ngit config --global github.user \"pnavaro\""
},

{
    "location": "slides/#Create-the-Julia-package-LyonCalcul-1",
    "page": "Talk",
    "title": "Create the Julia package LyonCalcul",
    "category": "section",
    "text": "(v1.0) pkg> generate LyonCalcul\nGenerating project LyonCalcul:\n    LyonCalcul/Project.toml\n    LyonCalcul/src/LyonCalcul.jlshell> cat LyonCalcul/Project.toml\nauthors = [\"Pierre Navaro <pierre.navaro@math.cnrs.fr>\"]\nname = \"LyonCalcul\"\nuuid = \"417a5b38-18da-11e9-35ce-9bdc85ad86c9\"\nversion = \"0.1.0\"\n\n[deps]shell> cat LyonCalcul/src/LyonCalcul.jl\nmodule LyonCalcul\n\ngreet() = print(\"Hello World!\")\n\nend # modulemodule LyonCalcul\n\nexport Advection\n\ninclude(\"bspline.jl\")\ninclude(\"advection.jl\")\n\nend # module"
},

{
    "location": "slides/#Activate-your-package-1",
    "page": "Talk",
    "title": "Activate your package",
    "category": "section",
    "text": "(v1.0) pkg> activate LyonCalcul\n\n(LyonCalcul) pkg> instantiate\n  Updating registry at `~/.julia/registries/General`\n  Updating git-repo `https://github.com/JuliaRegistries/General.git`\n Resolving package versions...\n  Updating `~/JuliaProjects/LyonCalcul/Project.toml`\n [no changes]"
},

{
    "location": "slides/#Add-dependencies-1",
    "page": "Talk",
    "title": "Add dependencies",
    "category": "section",
    "text": "(LyonCalcul) pkg> add FFTW\n Resolving package versions...\n  Updating `~/JuliaProjects/LyonCalcul/Project.toml`\n  [7a1cc6ca] + FFTW v0.2.4\n  Updating `~/JuliaProjects/LyonCalcul/Manifest.toml`\n  [621f4979] + AbstractFFTs v0.3.2\n  [b99e7846] + BinaryProvider v0.5.3\n  [34da2185] + Compat v1.4.0\n  [8f4d0f93] + Conda v1.1.1\n  ...\" The Manifest.toml allows someone to replicate the exact version of the dependencies that was recorded in the manifest on e.g. another machine. For a package that is to be used as a library, this is not useful.However, for an “application”, i.e. something at “top level” (say your julia code to do the simulations in a scientific paper) then it is likely useful to be able to replicate that exact state and the Manifest is thus useful to check in.\"Kristoffer Carlsson (Julia Computing)"
},

{
    "location": "slides/#Check-the-documentation-1",
    "page": "Talk",
    "title": "Check the documentation",
    "category": "section",
    "text": "\njulia> using LyonCalcul\n\nhelp?> Advection\nsearch: Advection\n\n  Advection(n, p, delta)\n\n  n :: Number of points.\n  p :: Spline degree.\n  delta :: space size step.\n\n  ───────────────────────────────────────────────────────────────────────────────────────\n\n  advection! = Advection( n, p, delta)\n  advection!( f, alpha )\n\n  Create a function to compute the interpolating spline\n  of degree p of odd\n  degree of a 1D function f on a periodic uniform mesh, at\n  all points x after a displacement alpha.\n  Input f type is Vector{Float64} and is updated inplace."
},

{
    "location": "slides/#Add-a-test-1",
    "page": "Talk",
    "title": "Add a test",
    "category": "section",
    "text": "cd LyonCalcul\nmkdir testadd file runtests.jlshell> cat test/runtests.jl\nusing Test\nusing LyonCalcul\n\n@testset \"Test advection sinus function\" begin\n\n     xmin, xmax, nx = 0.0, 2π, 128\n     dx = (xmax - xmin) / nx\n     x = range(xmin, stop=xmax, length=nx+1)[1:end-1] |> collect\n     f = sin.(x)\n     advection! = Advection( nx, 5, dx )\n     advection!( f, 0.5)\n     @test maximum(abs.(f .- sin.(x .- 0.5))) ≈ 0.0 atol = 1e-12\n\nend"
},

{
    "location": "slides/#Verify-the-test-1",
    "page": "Talk",
    "title": "Verify the test",
    "category": "section",
    "text": "(LyonCalcul) pkg> test\n   Testing LyonCalcul\n Resolving package versions...\nTest Summary:             | Pass  Total\nTest advection sinus function |    1      1\n   Testing LyonCalcul tests passed"
},

{
    "location": "slides/#Documentation-1",
    "page": "Talk",
    "title": "Documentation",
    "category": "section",
    "text": "https://github.com/JuliaDocs/Documenter.jljulia> using DocumenterTools\nshell> pwd\n/Users/navaro/JuliaProjects/LyonCalcul\njulia> DocumenterTools.generate(\"docs\")\n[ Info: name of package automatically determined to be `LyonCalcul`.\n[ Info: deploying documentation to `~/JuliaProjects/LyonCalcul/docs`\n[ Info: Generating .gitignore at /Users/navaro/JuliaProjects/LyonCalcul/docs/.gitignore\n[ Info: Generating make.jl at /Users/navaro/JuliaProjects/LyonCalcul/docs/make.jl\n[ Info: Generating Project.toml at /Users/navaro/JuliaProjects/LyonCalcul/docs/Project.toml\n[ Info: Generating src/index.md at /Users/navaro/JuliaProjects/LyonCalcul/docs/src/index.mdshell> cat docs/src/index.md\n# LyonCalcul.jl\n\nDocumentation for LyonCalcul.jl\n\n## Types and Functions```@autodocs\nModules = [LyonCalcul]\nOrder   = [:type, :function]\n```shell> cat docs/make.jl\n\nusing Documenter\nusing LyonCalcul\n\nmakedocs(modules=[LyonCalcul],\n         doctest = false,\n         format = Documenter.HTML(),\n         sitename = \"LyonCalcul.jl\",\n         pages = [\"Documentation\"    => \"index.md\"])\n\ndeploydocs(\n    deps   = Deps.pip(\"mkdocs\", \"python-markdown-math\"),\n    repo   = \"github.com/pnavaro/LyonCalcul.jl.git\",\n )"
},

{
    "location": "slides/#Add-a-repository-on-Github-1",
    "page": "Talk",
    "title": "Add a repository on Github",
    "category": "section",
    "text": "https://github.com/pnavaro/LyonCalcul.jlNote : the repository name has the \".jl\" extension$ echo \"# LyonCalcul.jl\" >> README.md\n$ git init\nInitialized empty Git repository in /Users/navaro/JuliaProjects/LyonCalcul/.git/\n$ git add .$ git commit -m \"first commit\"\n[master (root-commit) 8863c2e] first commit\n 11 files changed, 287 insertions(+)\n create mode 100644 Manifest.toml\n create mode 100644 Project.toml\n create mode 100644 README.md\n create mode 100644 docs/.gitignore\n create mode 100644 docs/Project.toml\n create mode 100644 docs/make.jl\n create mode 100644 docs/src/index.md\n create mode 100644 src/LyonCalcul.jl\n create mode 100644 src/advection.jl\n create mode 100644 src/bspline.jl\n create mode 100644 test/runtests.jl$ git remote add origin git@github.com:pnavaro/LyonCalcul.jl.git\n$ git push -u origin master\nEnumerating objects: 17, done.\nCounting objects: 100% (17/17), done.\nDelta compression using up to 8 threads\nCompressing objects: 100% (13/13), done.\nWriting objects: 100% (17/17), 4.29 KiB | 2.15 MiB/s, done.\nTotal 17 (delta 0), reused 0 (delta 0)\nTo github.com:pnavaro/LyonCalcul.jl.git\n * [new branch]      master -> master\nBranch \'master\' set up to track remote branch \'master\' from \'origin\'."
},

{
    "location": "slides/#Ignore-some-files-1",
    "page": "Talk",
    "title": "Ignore some files",
    "category": "section",
    "text": "$ cat .gitignore\n*.jl.cov\n*.jl.*.cov\n*.jl.mem\ndocs/build/\ndocs/site/\nManifest.toml"
},

{
    "location": "slides/#Install-first-version-of-Example-package-in-your-julia-installation-1",
    "page": "Talk",
    "title": "Install first version of Example package in your julia installation",
    "category": "section",
    "text": "(v1.1) pkg> add https://github.com/pnavaro/LyonCalcul.jl.git\n  Updating registry at `~/.julia/registries/General`\n  Updating git-repo `https://github.com/JuliaRegistries/General.git`\n   Cloning git-repo `https://github.com/pnavaro/LyonCalcul.jl.git`\n  Updating git-repo `https://github.com/pnavaro/LyonCalcul.jl.git`\n Resolving package versions...\n  Updating `~/.julia/environments/v1.1/Project.toml`\n  [417a5b38] ~ LyonCalcul v0.1.0 [`~/JuliaProjects/LyonCalcul`] ⇒ v0.1.0 #master (https://github.com/pnavaro/LyonCalcul.jl.git)\n  Updating `~/.julia/environments/v1.1/Manifest.toml`\n  [417a5b38] ~ LyonCalcul v0.1.0 [`~/JuliaProjects/LyonCalcul`] ⇒ v0.1.0 #master (https://github.com/pnavaro/LyonCalcul.jl.git`"
},

{
    "location": "slides/#Test-it-1",
    "page": "Talk",
    "title": "Test it",
    "category": "section",
    "text": "(v1.1) pkg> test LyonCalcul\n   Testing LyonCalcul\n\nTest Summary:                 | Pass  Total\nTest advection sinus function |    1      1\n   Testing LyonCalcul tests passed"
},

{
    "location": "slides/#Push-the-package-on-github-1",
    "page": "Talk",
    "title": "Push the package on github",
    "category": "section",
    "text": "cd LyonCalcul\necho \"# LyonCalcul.jl\" >> README.md\ngit init\ngit add README.md\ngit commit -m \"first commit\"\ngit remote add origin git@github.com:pnavaro/LyonCalcul.jl.git\ngit push -u origin master"
},

{
    "location": "slides/#On-Github-choose-your-license-1",
    "page": "Talk",
    "title": "On Github choose your license",
    "category": "section",
    "text": "Above the file list, click Create new file.In the file name field, type LICENSE (with all caps).Choose a license template button.\nClick Choose a license template.\nAdd a license to your project.\nDon\'t create pull request choose \"master\" branch."
},

{
    "location": "slides/#On-your-computer-1",
    "page": "Talk",
    "title": "On your computer",
    "category": "section",
    "text": "git pull origin master"
},

{
    "location": "slides/#Travis-1",
    "page": "Talk",
    "title": "Travis",
    "category": "section",
    "text": "Add your repository by going to https://travis-ci.orgProfile -> Settings -> Enable your repository"
},

{
    "location": "slides/#Codecov-1",
    "page": "Talk",
    "title": "Codecov",
    "category": "section",
    "text": "Add your repository by going to https://codecov.io/ghlanguage: julia\n\nos:\n  - linux\n  - osx\n\njulia:\n  - 1.0\n  - nightly\n\nnotifications:\n  email: true\n\nafter_success:\n    - julia -e \'using Pkg; cd(Pkg.dir(\"LyonCalcul\")); Pkg.add(\"Coverage\"); using Coverage; Codecov.submit(process_folder())\'\n\njobs:\n  include:\n    - stage: \"Documentation\"\n      julia: 1.0\n      os: linux\n      script:\n        - julia --project=docs/ -e \'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()\'\n        - julia --project=docs/ docs/make.jl\n      name: \"HTML\"\n      after_success: skip"
},

{
    "location": "slides/#Hosting-your-documentation-on-Github-pages-1",
    "page": "Talk",
    "title": "Hosting your documentation on Github pages",
    "category": "section",
    "text": "Launch julia in your package directory.pkg> add DocumenterTools\npkg> activate .julia> using DocumenterTools\njulia> using LyonCalcul\njulia> Travis.genkeys(LyonCalcul)Follow the instructions that are printed outAdd the public ssh key to your settings page for the GitHub repository.\nDon\'t forget to check Allow write access to allow Documenter to commit the   generated documentation to the repo.\nNext add the long private key named DOCUMENTER_KEY to the Travis settings page using the provided link.\nMake sure that the \"Display value in build log\" option is OFF."
},

{
    "location": "slides/#Update-documentation-with-Travis-1",
    "page": "Talk",
    "title": "Update documentation with Travis",
    "category": "section",
    "text": "To tell Travis that we want a new build stage we can add the following to the .travis.yml file:jobs:\n  include:\n    - stage: \"Documentation\"\n      julia: 1.0\n      os: linux\n      script:\n        - julia --project=docs/ -e \'using Pkg; Pkg.develop(PackageSpec(path=pwd()));\n                                               Pkg.instantiate()\'\n        - julia --project=docs/ docs/make.jl\n      after_success: skip"
},

{
    "location": "slides/#Enable-GitHub-Pages-1",
    "page": "Talk",
    "title": "Enable GitHub Pages",
    "category": "section",
    "text": "On GitHub, navigate to your GitHub Pages site\'s repository.\nUnder your repository name, click Settings.\nUse the Select source drop-down menu to select master or gh-pages as your GitHub Pages publishing source.\nClick Save.By default Documenter will create a link called dev that points to the latest versionhttps://pnavaro.github.io/LyonCalcul.jl/dev"
},

{
    "location": "slides/#Badges-1",
    "page": "Talk",
    "title": "Badges",
    "category": "section",
    "text": "It is common practice to make use of \"badges\" for Travis build status, code coverage and documentation. Adding the following to your package README.md should be all that is necessary:Codecov badge : https://codecov.io/gh/pnavaro/LyonCalcul.jl/settings/badge\nTravis badge : https://travis-ci.org/ click on the the badge corresponding to your repository.[![Build Status](https://travis-ci.org/pnavaro/LyonCalcul.jl.svg?branch=master)](https://travis-ci.org/pnavaro/LyonCalcul.jl)\n[![codecov](https://codecov.io/gh/pnavaro/LyonCalcul.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pnavaro/LyonCalcul.jl)\n[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pnavaro.github.io/LyonCalcul.jl/stable)\n[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pnavaro.github.io/LyonCalcul.jl/dev)"
},

{
    "location": "slides/#Register-your-package-1",
    "page": "Talk",
    "title": "Register your package",
    "category": "section",
    "text": "Set up AttoBot on your repository.\nYou need to tag your verson with git (for example v0.1.0)\nUse Github releases.\nWait a couple of days.I did not do it for LyonCalcul"
},

{
    "location": "slides/#Items-not-covered-1",
    "page": "Talk",
    "title": "Items not covered",
    "category": "section",
    "text": "Binary package PackageCompiler.jl\nMixed language BinDeps.jl\nJulia Observer\nCreate a pdf with Documenter\nLiterate.jl : create markdown file and/or jupyter notebook from a julia program. Easy way to create your examples and tutorials.\nWeaveAwayNotebooks.jl : convert Jupyter Notebooks to Literate.jl files\nWeave.jl : Scientific reports/literate programming for Julia"
},

{
    "location": "slides/#Bonus-1",
    "page": "Talk",
    "title": "Bonus",
    "category": "section",
    "text": "To set your documentation logo, just add a image file named logo.png in docs/src/assets directory.Its size must be 100 x 100 pixels.You can modify the julia logo images available on JuliaGraphics"
},

{
    "location": "slides/#Links-1",
    "page": "Talk",
    "title": "Links",
    "category": "section",
    "text": "Simplifying working with Julia packages and dependencies\nCreating a new package in Julia\nDocumenter\nFinalizing Your Julia Package: Documentation, Testing, Coverage, and Publishing (outdated)\nRevise.jl : Automatically update function definitions in a running Julia sessionThis page was generated using Literate.jl."
},

]}
