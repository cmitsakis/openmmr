# openMMR

Player *ranking* and *matchmaking* for team-based multiplayer games.

Player *ranking* is based on the [Weng-Lin Bayesian ranking algorithm](https://www.csie.ntu.edu.tw/~cjlin/papers/online_ranking/online_journal.pdf) which (unlike *Elo*) is suitable for multiplayer games.
The algorithm assumes the rating of each player is a random variable following a normal distribution with mean *μ*, and standard deviation *σ*.
Initially *μ* = 25 and *σ* = 25/3.
Each game updates those values so the mean (*μ*) converges to the actual strength of each player.
The actual strength of a player is within 2 standard deviations (*σ*) of the mean (*μ*) with probability 95%.
For example if *μ* = 29 and *σ* = 3 we can say we are 95% confident that the player's actual strength is between 23 and 35.

*Matchmaking* works by trying all possible teams of equal number of players.
The strength of a team is calculated by adding the ratings of all the players.

## Usage

Create a text file like this:

```
alice john mary > bob elisa thomas
alice bob mary < john elisa thomas
elisa john mary > bob thomas alice
elisa alice mary = bob thomas john
```

Each line represents a game. A > B means team A beat team B. A = B means the game was a draw.

Save the file as `games.txt`

#### Ranking

Run the following command to get ratings:

```sh
openmmr rank games.txt
```

output:

```
john    29.80  7.91
elisa   26.34  7.91
mary    26.14  7.91
thomas  23.86  7.91
alice   23.47  7.91
bob     20.39  7.91
```

The 2nd column is the mean (*μ*), and the 3rd column is the standard deviation (*σ*).

#### Matchmaking

Run the following command to get all the possible teams sorted from the most balanced to the most imbalanced.

```sh
openmmr teams games.txt
```

output:

```
[john thomas bob]  vs  [elisa mary alice]  74.05 vs 75.95  -1.90
[john mary bob]  vs  [elisa thomas alice]  76.32 vs 73.68  2.65
[john alice bob]  vs  [elisa mary thomas]  73.66 vs 76.34  -2.68
[john elisa bob]  vs  [mary thomas alice]  76.53 vs 73.47  3.06
[john thomas alice]  vs  [elisa mary bob]  77.13 vs 72.87  4.26
[john mary alice]  vs  [elisa thomas bob]  79.41 vs 70.59  8.81
[john elisa alice]  vs  [mary thomas bob]  79.61 vs 70.39  9.22
[john mary thomas]  vs  [elisa alice bob]  79.80 vs 70.20  9.59
[john elisa thomas]  vs  [mary alice bob]  80.00 vs 70.00  10.00
[john elisa mary]  vs  [thomas alice bob]  82.27 vs 67.73  14.55
```

You can also generate teams by using as input a file that contains ratings (instead of games) like this:

```sh
openmmr rank games.txt > ratings.txt
openmmr teams-from-rank ratings.txt
```

Or you can use "-" to read from the standard input like this:
```sh
openmmr rank games.txt | openmmr teams-from-rank -
```

## Supported platforms

*Linux*, *Windows 8.1+*, *macOS 10.13+*

## Installation

### Option 1: Download release binary (recommended)

Download the latest [release](https://github.com/cmitsakis/openmmr/releases) and run it. No installation is required.

#### macOS

On *macOS* you have to remove the application from *quarantine* by following the instructions [here](https://support.apple.com/guide/mac-help/welcome/mac), or by running the following command:

```sh
xattr -d com.apple.quarantine /path/to/openmmr
```

### Option 2: Build from source

Requires [Zig](https://ziglang.org/) v0.8.1

Run `zig build` at the root of the repository.

## License

Copyright (C) 2021 Charalampos Mitsakis

*openMMR* is licensed under the [EUPL-1.2-or-later](LICENSE).

The implementation of the *Weng-Lin Bayesian ranking algorithm* in the file `weng_lin/rank.c` is licensed under the terms of the [BSD-3-Clause](weng_lin/LICENSE) license.
