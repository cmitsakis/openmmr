// Copyright (C) 2021 Charalampos Mitsakis
// Licensed under the EUPL-1.2-or-later

const std = @import("std");
const log = std.log;
const mem = std.mem;
const rating = @import("rating.zig");
const matchmaking = @import("matchmaking.zig");

const usage =
    \\USAGE: openmmr <subcommand> [<file>]
    \\
    \\SUBCOMMANDS:
    \\  rank:            calculate rating for every player
    \\  teams:           generate all the possible teams sorted from the most balanced to the most imbalanced
    \\  teams-from-rank: generate teams from rank
    \\  help:            print help
    \\  license:         print license
;

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = &arena.allocator;

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

    if (args.len != 3) {
        try stderr.print("error: too many arguments\n", .{});
        try stderr.print("{s}\n", .{usage});
        return;
    }

    const filename = args[2];
    const file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();
    const reader = file.reader();
    if (mem.eql(u8, cmd, "rank")) {
        const score_by_player = try rating.scoreByPlayerFromGames(allocator, reader);
        const players = try rating.scoreByPlayerToSliceOfPlayers(allocator, score_by_player);
        try rating.printRatings(players);
    } else if (mem.eql(u8, cmd, "teams-from-rank")) {
        const players = try matchmaking.readScores(allocator, reader);
        if (players.len > 64) {
            try stderr.print("error: too many players\n", .{});
            return;
        }
        try matchmaking.teams(players);
    } else if (mem.eql(u8, cmd, "teams")) {
        const score_by_player = try rating.scoreByPlayerFromGames(allocator, reader);
        const players = try rating.scoreByPlayerToSliceOfPlayers(allocator, score_by_player);
        if (players.len > 64) {
            try stderr.print("error: too many players\n", .{});
            return;
        }
        try matchmaking.teams(players);
    } else {
        try stderr.print("error: unknown command", .{});
        try stderr.print("{s}\n", .{usage});
    }
}
