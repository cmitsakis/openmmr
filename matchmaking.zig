// Copyright (C) 2021 Charalampos Mitsakis
// Licensed under the EUPL-1.2-or-later

const std = @import("std");
const log = std.log;
const mem = std.mem;
const types = @import("types.zig");

const Assignment = struct {
    binary: u64,
    team_0_rating_mu: f64,
    team_1_rating_mu: f64,
    fn lessThan(_: void, lhs: Assignment, rhs: Assignment) bool {
        return std.math.absFloat(lhs.team_0_rating_mu - lhs.team_1_rating_mu) < std.math.absFloat(rhs.team_0_rating_mu - rhs.team_1_rating_mu);
    }
    fn ratingDiff(self: Assignment) f64 {
        return self.team_0_rating_mu - self.team_1_rating_mu;
    }
};

pub fn printTeams(players: []types.Player) !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = &arena.allocator;

    if ((players.len / 2) * 2 != players.len) {
        return error.OddNumOfPlayers;
    }

    var total_rating_mu: f64 = 0;
    for (players) |player| {
        total_rating_mu += player.rating.mu;
    }

    const u64_max: u64 = 0b1111111111111111111111111111111111111111111111111111111111111111;
    const b_start = ~(u64_max << @intCast(u6, players.len / 2));
    const b_end_mask = ~(u64_max << @intCast(u6, players.len));
    const b_end = ~b_start & b_end_mask;
    //log.debug("b_start: {b} b_end: {b}", .{ b_start, b_end });
    var team_assignments = std.ArrayList(Assignment).init(allocator);
    var buf: [64]u8 = undefined;
    var b = b_start;
    while (b <= b_end) : (b = nextBinary(b)) {
        if (b & 1 == 1) {
            continue;
        }
        //log.debug("b: {b}", .{b});
        var team_1_rating_mu: f64 = 0;
        buf = undefined;
        var bString = try std.fmt.bufPrint(buf[0..], "{b}", .{b});
        for (buf) |bChar, i| {
            if (bChar == '1') {
                team_1_rating_mu += players[bString.len - 1 - i].rating.mu;
            }
        }
        var team_0_rating_mu: f64 = total_rating_mu - team_1_rating_mu;
        try team_assignments.append(Assignment{ .team_0_rating_mu = team_0_rating_mu, .team_1_rating_mu = team_1_rating_mu, .binary = b });
    }
    var team_assignments_slice = team_assignments.toOwnedSlice();
    std.sort.sort(Assignment, team_assignments_slice, {}, Assignment.lessThan);
    const stdout = std.io.getStdOut().writer();
    for (team_assignments_slice) |assignment| {
        var arena_assignment = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena_assignment.deinit();
        const allocator_assignment = &arena_assignment.allocator;
        var team_0_players = std.ArrayList(types.Username).init(allocator_assignment);
        var team_1_players = std.ArrayList(types.Username).init(allocator_assignment);
        b = assignment.binary;
        var i: usize = 0;
        while (i < players.len) {
            if (b & 1 == 0) {
                try team_0_players.append(players[i].username);
            } else {
                try team_1_players.append(players[i].username);
            }
            b = b >> 1;
            i = i + 1;
        }
        try stdout.print("[{s}]", .{std.mem.join(allocator_assignment, " ", team_0_players.toOwnedSlice())});
        try stdout.print("  vs  ", .{});
        try stdout.print("[{s}]", .{std.mem.join(allocator_assignment, " ", team_1_players.toOwnedSlice())});
        try stdout.print("  {d:.2} vs {d:.2}  {d:.2}\n", .{ assignment.team_0_rating_mu, assignment.team_1_rating_mu, assignment.ratingDiff() });
    }
}

pub fn readScores(allocator: *std.mem.Allocator, reader: std.fs.File.Reader) ![]types.Player {
    var players = std.ArrayList(types.Player).init(allocator);
    var line_buf: [1024]u8 = undefined;
    while (reader.readUntilDelimiterOrEof(&line_buf, '\n') catch |err| {
        return err;
    }) |line| {
        var tokens = std.mem.tokenize(line, " ");
        const player_username = try std.mem.dupe(allocator, u8, tokens.next().?);
        const player_rating_mu = try std.fmt.parseFloat(f64, tokens.next().?);
        const player_rating_sigma = try std.fmt.parseFloat(f64, tokens.next().?);
        const p: types.Player = types.Player{
            .username = player_username,
            .rating = types.Rating{
                .mu = player_rating_mu,
                .sigma = player_rating_sigma,
            },
        };
        try players.append(p);
    }
    return players.toOwnedSlice();
}

fn nextBinary(v: u64) u64 {
    // copied from: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    const t: u64 = (v | (v -% 1)) +% 1;
    return t | ((((t & (0 -% t)) / (v & (0 -% v))) >> 1) - 1);
}

const expect = @import("std").testing.expect;

test "nextBinary" {
    var n: u64 = 0b00010011;
    n = nextBinary(n);
    try expect(n == 0b00010101);
    n = nextBinary(n);
    try expect(n == 0b00010110);
    n = nextBinary(n);
    try expect(n == 0b00011001);
    n = nextBinary(n);
    try expect(n == 0b00011010);
    n = nextBinary(n);
    try expect(n == 0b00011100);
    n = nextBinary(n);
    try expect(n == 0b00100011);

    n = 0b01110010;
    n = nextBinary(n);
    try expect(n == 0b01110100);
    n = nextBinary(n);
    try expect(n == 0b01111000);
    n = nextBinary(n);
    try expect(n == 0b10000111);
    n = nextBinary(n);
    try expect(n == 0b10001011);
    n = nextBinary(n);
    try expect(n == 0b10001101);
    n = nextBinary(n);
    try expect(n == 0b10001110);
}
