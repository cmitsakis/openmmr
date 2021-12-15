// Copyright (C) 2021 Charalampos Mitsakis
// Licensed under the EUPL-1.2-or-later

const std = @import("std");
const log = std.log;
const mem = std.mem;
const types = @import("types.zig");
const wengLin = @cImport(@cInclude("weng_lin/rank.h"));

const init_mu: f64 = 25;
const init_sigma: f64 = init_mu / 3;
const beta: f64 = init_sigma * 0.5;
const num_teams: c_int = 2;

const Game = struct {
    players: [2][]types.Username,
    rank: [2]c_int,
};

const Algorithm = enum {
    BT_full,
    //BT_partial,
    PL,
    TM_full,
};

pub fn scoreByPlayerFromGames(allocator: *std.mem.Allocator, reader: std.fs.File.Reader) !std.StringHashMap(types.Rating) {
    var score_by_player = std.StringHashMap(types.Rating).init(allocator);

    var line_buf: [1024]u8 = undefined;
    while (reader.readUntilDelimiterOrEof(&line_buf, '\n') catch |err| {
        return err;
    }) |line| {
        if (line[0] == '#') {
            continue;
        }
        if (line.len == 0) {
            continue;
        }
        var arena_line = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena_line.deinit();
        const allocator_line = &arena_line.allocator;
        const game = try parseGame(allocator_line, line);
        for (game.players) |team_usernames| {
            for (team_usernames) |username| {
                if (!score_by_player.contains(username)) {
                    const username2 = try std.mem.dupe(allocator, u8, username); // use 'allocator' instead of 'allocator_line' because usernames outlive this function
                    try score_by_player.put(username2, types.Rating{ .mu = init_mu, .sigma = init_sigma });
                }
            }
        }
        try processGame(allocator_line, &score_by_player, game, Algorithm.BT_full);
    }
    return score_by_player;
}

pub fn scoreByPlayerToSliceOfPlayers(allocator: *std.mem.Allocator, score_by_player: std.StringHashMap(types.Rating)) ![]types.Player {
    var players_rated = std.ArrayList(types.Player).init(allocator);
    var score_by_player_iterator = score_by_player.iterator();
    while (score_by_player_iterator.next()) |e| {
        try players_rated.append(types.Player{ .username = e.key_ptr.*, .rating = e.value_ptr.* });
    }
    var players_rated_slice = players_rated.toOwnedSlice();
    std.sort.sort(types.Player, players_rated_slice, {}, types.Player.ratingGreaterThan);
    return players_rated_slice;
}

fn parseGame(allocator: *std.mem.Allocator, line: []const u8) !Game {
    var team_0 = std.ArrayList(types.Username).init(allocator);
    var team_1 = std.ArrayList(types.Username).init(allocator);
    var tokens = std.mem.tokenize(line, " ");
    var team_seperator_seen = false;
    var rank: [2]c_int = undefined;
    while (tokens.next()) |token| {
        if (std.mem.eql(u8, token, ">") or std.mem.eql(u8, token, "<") or std.mem.eql(u8, token, "=")) {
            if (std.mem.eql(u8, token, ">")) {
                rank = [_]c_int{ 0, 1 };
            } else if (std.mem.eql(u8, token, "<")) {
                rank = [_]c_int{ 1, 0 };
            } else if (std.mem.eql(u8, token, "=")) {
                rank = [_]c_int{ 0, 0 };
            }
            if (!team_seperator_seen) {
                team_seperator_seen = true;
            } else {
                return error.InvalidChar;
            }
            continue;
        } else {
            if (team_seperator_seen) {
                try team_1.append(token);
            } else {
                try team_0.append(token);
            }
        }
    }
    if (!team_seperator_seen) {
        return error.InvalidChar;
    }
    return Game{
        .players = [2][]types.Username{ team_0.toOwnedSlice(), team_1.toOwnedSlice() },
        .rank = rank,
    };
}

fn processGame(allocator: *std.mem.Allocator, score_by_player: *std.StringHashMap(types.Rating), game: Game, algo: Algorithm) !void {
    var rank: [2]c_int = game.rank; // copy rank to a mutable var in order to use @ptrCast() later
    var teams_players_usernames: [2][]types.Username = undefined;
    var teams_players_mu: [2][*c]f64 = undefined;
    var teams_players_sigma: [2][*c]f64 = undefined;
    var teams_num_players: [2]c_int = undefined;
    for (game.players) |team_usernames, i| {
        var team_players_mu = std.ArrayList(f64).init(allocator);
        var team_players_sigma = std.ArrayList(f64).init(allocator);
        var team_players_usernames = std.ArrayList(types.Username).init(allocator);
        for (team_usernames) |username| {
            try team_players_usernames.append(username);
            const val = score_by_player.get(username).?;
            try team_players_mu.append(val.mu);
            try team_players_sigma.append(val.sigma);
        }
        var team_players_mu_slice = team_players_mu.toOwnedSlice();
        var team_players_sigma_slice = team_players_sigma.toOwnedSlice();
        teams_players_usernames[i] = team_players_usernames.toOwnedSlice();
        teams_players_mu[i] = @ptrCast([*c]f64, team_players_mu_slice.ptr);
        teams_players_sigma[i] = @ptrCast([*c]f64, team_players_sigma_slice.ptr);
        teams_num_players[i] = @intCast(c_int, team_players_mu_slice.len);
    }
    var predictions: [30]c_int = undefined;
    _ = switch (algo) {
        .BT_full => wengLin.update_BT_full(beta, num_teams, @ptrCast([*c]c_int, &rank), @ptrCast([*c]c_int, &teams_num_players), @ptrCast([*c][*c]f64, &teams_players_mu), @ptrCast([*c][*c]f64, &teams_players_sigma), @ptrCast([*c]c_int, &predictions)),
        .PL => wengLin.update_PL(beta, num_teams, @ptrCast([*c]c_int, &rank), @ptrCast([*c]c_int, &teams_num_players), @ptrCast([*c][*c]f64, &teams_players_mu), @ptrCast([*c][*c]f64, &teams_players_sigma), @ptrCast([*c]c_int, &predictions)),
        .TM_full => wengLin.update_TM_full(beta, num_teams, @ptrCast([*c]c_int, &rank), @ptrCast([*c]c_int, &teams_num_players), @ptrCast([*c][*c]f64, &teams_players_mu), @ptrCast([*c][*c]f64, &teams_players_sigma), @ptrCast([*c]c_int, &predictions)),
    };
    for (teams_players_usernames) |team_players_usernames, i| {
        for (team_players_usernames) |username, j| {
            var val = score_by_player.get(username).?;
            val.mu = teams_players_mu[i][j];
            val.sigma = teams_players_sigma[i][j];
            try score_by_player.put(username, val);
        }
    }
}

pub fn printRatings(players_rated_slice: []types.Player) !void {
    comptime var max_padding_string: []const u8 = "                                ";
    comptime var max_padding_len: usize = 32;
    const stdout = std.io.getStdOut().writer();
    // find max username width
    var max_width_username: usize = 0;
    for (players_rated_slice) |player| {
        if (player.username.len > max_width_username) {
            max_width_username = player.username.len;
        }
    }
    // find max rating width
    var buf: [max_padding_len]u8 = undefined;
    var max_width_rating: usize = 0;
    for (players_rated_slice) |player| {
        const l = lengthFloat(buf[0..], player.rating.mu);
        if (l > max_width_rating) {
            max_width_rating = l;
        }
    }
    // find max sigma width
    var max_width_sigma: usize = 0;
    for (players_rated_slice) |player| {
        const l = lengthFloat(buf[0..], player.rating.sigma);
        if (l > max_width_sigma) {
            max_width_sigma = l;
        }
    }
    // print for each player
    for (players_rated_slice) |player| {
        const padding_username = reduceBelow(max_padding_len, max_width_username - player.username.len);
        const padding_rating = reduceBelow(max_padding_len, max_width_rating - lengthFloat(buf[0..], player.rating.mu));
        const padding_sigma = reduceBelow(max_padding_len, max_width_sigma - lengthFloat(buf[0..], player.rating.sigma));
        try stdout.print("{s}{s}  {s}{d:>.2}  {s}{d:.2}\n", .{ player.username, max_padding_string[0..padding_username], max_padding_string[0..padding_rating], player.rating.mu, max_padding_string[0..padding_sigma], player.rating.sigma });
    }
}

fn reduceBelow(max: usize, s: usize) usize {
    if (s < max) {
        return s;
    } else {
        return max - 1;
    }
}

fn lengthFloat(buf: []u8, f: f64) usize {
    if (std.fmt.bufPrint(buf, "{d:.2}", .{f})) |str| {
        return str.len;
    } else |err| switch (err) {
        error.NoSpaceLeft => {
            return buf.len;
        },
    }
    return str.len;
}
