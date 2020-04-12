## How to contribute to SpectralWaveData

### **Did you find a bug?**

* Please, ensure the bug was not already reported by searching this repository under 
  [Issues](https://github.com/SpectralWaveData/spectral_wave_data/issues).

* If you're unable to find an open issue addressing the problem, 
  [open a new one](https://github.com/SpectralWaveData/spectral_wave_data/issues/new). 
  Be sure to include a title, clear description, and as much 
  relevant information as possible. 
  When relevant please include code or an executable test case 
  demonstrating the expected behavior that is not occurring.

### **Did you write a patch that fixes a bug?**

* For any substantial change of code, you should first make a
  [fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) of this repository and implement the related changes 
  in a new branch of this fork.
* When you have updated and tested your new code, open a new GitHub
  pull request. Follow this GitHub 
  [guideline](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests)
  on how to collaborate with issues and pull requests.
  
* Ensure the pull request description clearly describes the problem and
  solution. 
  Include the relevant issue number if applicable.

### **When will my contribution be included in the official version?**

It depends on how extensive and critical your contribution is, and on 
available resources. In some cases, we may also have to adapt your contribution. 
Only acceptable contributions will be included. Please, discuss in the
community to check if the contribution is acceptable before you start
to implement some major new contribution. 

This is the general procedure we follow if we receive a pull request from
you. We basically follow the 
[gitflow](https://datasift.github.io/gitflow/IntroducingGitFlow).

* If your contribution is related to the documentation only,
  and it is important for clarification of the current version, like fixing
  a typo in an equation, 
  we will first update the 'hotfix-docs' branch. 
  After some testing, this branch will be merged into the 'master' branch.
  No new version number and libraries will be created, 
  but the official repository and the 'latest'-version of the live 
  documentation will be updated. A note will be inserted in the release log.

* For change of code, 
  we will check out your branch and double check that you have enough 
  test code and documentation for your new feature. 
  If some of your test code, or existing test code does not work, or
  the new feature lack documentation, we will not proceed. We let you know!

* For minor but urgent fixes we will merge your contribution into a new hotfix
  branch, soon to be merged into the official 'master' branch after positive
  testing. Then a new official release, including compiled binaries, 
  will be announced in the 'release' tab of this repository. Make sure
  you have set notification for release announcements.

* For other contributions, we will try to merge your branch into the
  'develop' branch in our official repository. This branch may
  include other 'finalized' features too. The software version tag 
  in this branch will include an alpha tag. We will let you know when
  your feature has been successfully merged into this branch. We would
  be happy if you review the related code and test it to see if it solves
  your original purpose.

* When the 'develop' branch becomes mature, when all the pytest 
  code runs successfully and the related documentation is properly
  integrated, a new feature complete 'release' branch is established. 
  The new version number in this branch will include a beta tag. We 
  will let you know when your code is included in a release branch.

* Only bug fixes and clarifications of the documentation,
  will be made in the release branch. It is highly 
  appreciated if you test the release version on your production cases
  and let us know your findings. Negative or positive.

* When the release branch has matured, we will merge it into the 'master'
  branch for a new official release. When all pytests are ok, we will
  push the release to GitHub, upload new binaries, and announce it 
  in the 'release' tab of this repository.

* If major updates are introduced in a new version, a release candidate
  may first be released in the 'master' branch. 
  E.g. 1.1.0-rc1. Typical a month later 1.1.0.

The reason for this approach is that for many companies, the quality
of the 'master' branch is highly critical for their production. 
For more research-oriented users, the 'feature', 'develop' and 
'release' branches may be more interesting.



### **There are no shape classes supporting the formulation in my wave generator. What should I do?**

A main goal with this API is to facilitate utilization of all spectral
formulations - including yours.

To obtain that goal, in addition to minimize the maintenance 
cost and the cost for vendors to develop highly optimized SWD
implementations, it is important to apply generic shape classes that can
be reused by many wave generators. Consequently, _if_ it is possible to 
reuse an existing shape class for your formulation, we will most likely
not include an alternative class.

Mathematics is by nature a generic concept. There are not many
fundamental spectral solutions of Laplace's equation.
As an example, assuming the spacing of wave numbers is constant, any long-crested wave fields periodic in the horizontal 
space, can be expressed using shape class 3 as defined in the 
theory documentation. This is independent on 
your free surface and sea floor boundary conditions. Infinite depth or highly 
varying bathymetry. You may apply specific hyperbolic functions or splitting of the potential in your wave generator, but the output of your
total wave field from your wave generator can always be adjusted to fit 
the sum of exponential shape functions as defined in shape-3. If not, you are not solving Laplace's equation. Explicit transformations are discussed in 
the documentation, in the theory section for shape-3. 

This does not mean that you should reformulate your wave generator. 
Just keep your formulation and implementation as it is, but in your new output
routine you should apply minor transformations (if needed) when writing
to the SWD file. This is a small price to pay, compared to maintaining many 
more or less identical shape classes.
The documentation presents two case studies on how to relate
specific real-world wave generators to the generic formulation for writing
SWD compliant output. (See section 7.1.1)

Yes, there are sound formulations that cannot be represented by the
current set of shape classes. One such example is short crested waves
propagating over a varying seabed. We would be happy if you would
like to implement such a class or in some other way help
us on developing such a class. Ideas, testing etc. Before starting to 
write theory documentation, implement the code and test suites etc.,
the community should first discuss the formulation to ensure we implement
an efficient generic formulation suitable for many wave generators.

If you still are unsure how to provide SWD output from your wave generator
we may provide further help. If you send us a document defining your 
coordinate system, time domain velocity potential and surface elevation, 
we may suggest specific transformations and shape classes for your
considerations. We look forward to seeing your wave generator producing
SWD files, ready for a growing set of applications.


### **Do you have other suggestions?**

Interesting, please let us know...


Best regards,

The SpectralWaveData Team
