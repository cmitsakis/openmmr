// Copyright (C) 2021 Charalampos Mitsakis
// Licensed under the EUPL-1.2-or-later

pub const Username = []const u8;

pub const Player = struct {
    username: Username,
    rating: Rating,
    pub fn ratingGreaterThan(_: void, lhs: Player, rhs: Player) bool {
        return lhs.rating.mu > rhs.rating.mu;
    }
};

pub const Rating = struct {
    mu: f64,
    sigma: f64,
};
