const Builder = @import("std").build.Builder;

const cflags = [_][]const u8{
    "-std=c99",
    "-pedantic",
    "-Werror",
    "-Wall",
};

pub fn build(b: *Builder) void {
    const exe = b.addExecutable("openmmr", "main.zig");
    exe.setTarget(b.standardTargetOptions(.{}));
    exe.setBuildMode(b.standardReleaseOptions());
    exe.addIncludeDir(".");
    exe.addCSourceFile("weng_lin/rank.c", &cflags);
    exe.linkSystemLibrary("c");

    b.default_step.dependOn(&exe.step);
    b.installArtifact(exe);
}
