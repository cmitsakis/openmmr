// Copyright (C) 2021 Charalampos Mitsakis
// Licensed under the EUPL-1.2-or-later

const std = @import("std");
const log = std.log;
const mem = std.mem;
const ranking = @import("ranking.zig");
const matchmaking = @import("matchmaking.zig");
const types = @import("types.zig");

const usage =
    \\USAGE: openmmr <subcommand> [<file>]
    \\
    \\SUBCOMMANDS:
    \\  rank:            calculate ratings for every player
    \\  teams:           generate all the possible teams sorted from the most balanced to the most imbalanced
    \\  teams-from-rank: generate teams from rank
    \\  help:            print help
    \\  license:         print license
;

const emptyStruct = struct {};

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const args = try std.process.argsAlloc(allocator);

    const stdout = std.io.getStdOut().writer();
    const stderr = std.io.getStdErr().writer();

    if (args.len <= 1) {
        try stderr.print("error: no arguments\n", .{});
        try stderr.print("{s}\n", .{usage});
        return;
    }

    const cmd = args[1];

    if (mem.eql(u8, cmd, "license")) {
        try stdout.print("Copyright (C) 2021 Charalampos Mitsakis\n", .{});
        try stdout.print("openMMR is licensed under the EUPL-1.2-or-later\n", .{});
        try stdout.print("\n{s}\n", .{@embedFile("LICENSE")});
        try stdout.print("\n[3rd party licenses]\n\n", .{});
        try stdout.print("The implementation of the Weng-Lin Bayesian ranking algorithm\n", .{});
        try stdout.print("is licensed under the terms of the BSD-3-Clause license.\n", .{});
        try stdout.print("\n{s}\n", .{@embedFile("weng_lin/LICENSE")});
        return;
    }
    if (mem.eql(u8, cmd, "help")) {
        try stderr.print("{s}\n", .{usage});
        return;
    }

    var players_included = std.StringHashMap(emptyStruct).init(allocator);
    var players_excluded = std.StringHashMap(emptyStruct).init(allocator);
    if (args.len > 3) {
        var included_or_excluded: i32 = 0;
        var i: usize = 2;
        while (i < args.len - 1) : (i += 1) {
            const arg = args[i];
            if (mem.eql(u8, arg, "--include") or mem.eql(u8, arg, "--exclude")) {
                if (included_or_excluded != 0) {
                    try stderr.print("cannot use both --include and --exclude options\n", .{});
                    return;
                }
                included_or_excluded = if (mem.eql(u8, arg, "--include")) 1 else -1;
                i += 1;
                if (i == args.len - 1) {
                    try stderr.print("no file argument given\n", .{});
                    return;
                }
                var usernames = std.mem.tokenize(u8, args[i], ",");
                while (usernames.next()) |username| {
                    switch (included_or_excluded) {
                        1 => try players_included.put(username, .{}),
                        -1 => try players_excluded.put(username, .{}),
                        else => {
                            unreachable;
                        },
                    }
                }
            } else {
                try stderr.print("invalid argument: {s}\n", .{args[i]});
                return;
            }
        }
    }

    const filename = args[args.len - 1];
    // open file or stdin (if filename == "-")
    const file = if (!mem.eql(u8, filename, "-")) try std.fs.cwd().openFile(filename, .{}) else std.io.getStdIn();
    defer if (!mem.eql(u8, filename, "-")) file.close();
    const reader = file.reader();

    if (mem.eql(u8, cmd, "rank")) {
        const score_by_player = try ranking.scoreByPlayerFromGames(allocator, reader);
        const players_unfiltered = try ranking.scoreByPlayerToSliceOfPlayers(allocator, score_by_player);
        const players = try filterPlayersExclude(allocator, try filterPlayersInclude(allocator, players_unfiltered, players_included), players_excluded);
        try ranking.printRatings(players);
    } else if (mem.eql(u8, cmd, "teams-from-rank")) {
        const players_unfiltered = try matchmaking.readScores(allocator, reader);
        const players = try filterPlayersExclude(allocator, try filterPlayersInclude(allocator, players_unfiltered, players_included), players_excluded);
        if (players.len > 64) {
            try stderr.print("error: too many players\n", .{});
            return;
        }
        try matchmaking.printTeams(players);
    } else if (mem.eql(u8, cmd, "teams")) {
        var players_included_iterator = players_included.iterator();
        while (players_included_iterator.next()) |p| {
            try stdout.print("included: {s}\n", .{p.key_ptr.*});
        }
        const score_by_player = try ranking.scoreByPlayerFromGames(allocator, reader);
        const players_unfiltered = try ranking.scoreByPlayerToSliceOfPlayers(allocator, score_by_player);
        const players = try filterPlayersExclude(allocator, try filterPlayersInclude(allocator, players_unfiltered, players_included), players_excluded);
        if (players.len > 64) {
            try stderr.print("error: too many players\n", .{});
            return;
        }
        matchmaking.printTeams(players) catch |err| switch (err) {
            error.OddNumOfPlayers => try stderr.print(
                \\error: input file contains an odd number of players.
                \\Use the --exclude or --include options to get an even number of players like this:
                \\openmmr teams --exclude john,elisa games.txt
                \\openmmr teams --include alice,bob,thomas,mary games.txt
                \\
            , .{}),
            else => return err,
        };
    } else {
        try stderr.print("error: unknown command", .{});
        try stderr.print("{s}\n", .{usage});
    }
}

fn filterPlayersInclude(allocator: std.mem.Allocator, players: []types.Player, included: std.StringHashMap(emptyStruct)) ![]types.Player {
    if (included.count() == 0) {
        return players;
    }
    var players2 = std.ArrayList(types.Player).init(allocator);
    for (players) |player| {
        if (included.get(player.username)) |_| {
            try players2.append(player);
        }
    }
    return players2.toOwnedSlice();
}

fn filterPlayersExclude(allocator: std.mem.Allocator, players: []types.Player, excluded: std.StringHashMap(emptyStruct)) ![]types.Player {
    if (excluded.count() == 0) {
        return players;
    }
    var players2 = std.ArrayList(types.Player).init(allocator);
    for (players) |player| {
        if (excluded.get(player.username)) |_| {} else {
            try players2.append(player);
        }
    }
    return players2.toOwnedSlice();
}
