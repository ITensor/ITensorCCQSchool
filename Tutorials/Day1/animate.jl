"""
    animate(f; nframes::Int, fps::Int = 30)

Call `f(i)` for each frame `i in 1:nframes`, where `f(i)` returns an object to print to the
terminal at each frame. Renders each frame in-place in the terminal at `fps` frames per
second.
"""
function animate(f; nframes::Int, fps::Real = 30)
    io = IOBuffer()
    # Hide cursor
    print(stdout, "\x1b[?25l")
    return try
        # Clear screen
        print(stdout, "\x1b[2J")
        for i in 1:nframes
            seekstart(io)
            show(io, MIME("text/plain"), f(i))
            frame_str = String(take!(io))
            # Move home
            print(stdout, "\x1b[H")
            println(frame_str)
            sleep(1 / fps)
        end
    finally
        # Show cursor again
        print(stdout, "\x1b[?25h")
    end
end
