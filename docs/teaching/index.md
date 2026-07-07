# Teaching

This area collects student-facing teaching material built around OpenShieldHIT.
It is separate from the user manual and from runnable examples:

- user documentation explains the input format and command line;
- examples show complete reusable simulations;
- teaching material gives exercises, prompts, interpretation tasks, and lecture scaffolding.

The material is intended for graduate-level courses and workshops in Monte Carlo
transport, particle therapy physics, and related medical-physics topics.  The
first course is practical and tool-oriented; later courses can stand alone as
more focused physics or mathematics modules.

OpenShieldHIT is conceptually inspired by the same application domain and input
style as SHIELD-HIT12A, but it shares no source code with SHIELD-HIT or
SHIELD-HIT12A.

!!! note "Teaching-page formatting"
    Use admonition boxes sparingly and consistently:

    - `note` for assumptions and context;
    - `tip` for practical workflow hints;
    - `warning` for current limitations or common mistakes;
    - `abstract` for deliverables, project contracts, or grading criteria.

## Courses

| Course | Focus |
|---|---|
| [00 Getting started](courses/00-getting-started/index.md) | First simulations, input files, scoring, and plotting |
| [01 Basic Monte Carlo transport](courses/01-basic-monte-carlo-transport/index.md) | Histories, random numbers, sampling, estimators, and uncertainty |
| [02 Beam dynamics](courses/02-beam-dynamics/index.md) | Phase space, beam envelopes, spot lists, and Twiss parameters |
| [03 Atomic interactions](courses/03-atomic-interactions/index.md) | Stopping power, range, straggling, and multiple scattering |
| [04 Nuclear interactions](courses/04-nuclear-interactions/index.md) | Nuclear reactions, fragments, neutrons, and model comparison |

## Group projects

Larger problem-based-learning assignments live in
[projects/](projects/index.md): multi-week group projects that — unlike the
exercise sets — modify the code itself (new physics channels, new nuclear
benchmarks, performance and correctness engineering, detector-response
scorers).  They are written for graduate and post-graduate groups working
together with LLM coding agents, with scope tiers, validation requirements,
and stable implementation themes.

## Slides

Shared slide placeholders live in [slides/](slides/index.md).  Course-specific
slides live below each course directory.

## Exercise style

Exercise sets should be written so a student can work through them without
editing source code.  Prefer small, reproducible case directories and short
analysis scripts over manual steps that are hard to check.

Each exercise set should eventually contain:

- learning goals;
- required background;
- input files or links to prepared cases;
- tasks;
- expected observations;
- optional extension questions.
